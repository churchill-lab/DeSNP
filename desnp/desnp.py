# -*- coding: utf-8 -*-

import sys
import logging
import time
import csv
import gzip
from builtins import range

import pysam
from datetime import datetime
from desnp.probe import parse_probe


"""
desnp.py
January 9, 2012
Dave Walton - dave.walton@jax.org

August 12, 2019
Matt Vincent - matt.vincent@jax.org

This program is part of a toolkit for manipulating expression data.
This specific piece of code will take a file of probes, a set of strains,
and a file of SNPS and filter out all the probes that have a SNP in
one or more of the given strains against a reference sequence.  The
C57BL6/J from NCBI Build 37 (mm9) is assumed to be the reference.

Output will be written to standard out, and will be a file identical to the
probes input file (either <probes.txt> below or the 'probes.tsv' file if using
the -z option), except that the probes with snps will be either filtered out
of the file, or (if the -k option is used below) the probes will still be in the
file with a column added stating 'SNPS', in the case where snps are found in the
probe.

Usage:
  ./desnp.py [OPTIONS] -f <probes.txt> -g <snps.gz> -s <strains> (1st form)
  ./desnp.py [OPTIONS] -z <probes.zip> -g <snps.gz> -s <strains> (2nd form)
 
OPTIONS:
    -c, --comma    the probe file is comma delimited
    -f, --file     the probe file. This or -z are required
    -g, --gzipsnp  the gzipped snp file.  This requires an associated .tbi tabix
                   index file to be present in the same location
    -h, --help     return this message
    -i, --idcol    the name of the unique probe id column.  if not provided assumes 'id'
    -l, --log      same as verbose but sends diagnostics to desnp.log
    -o, --out      the name of the output file the results will go to
    -r, --returnstrains Can be used in conjunction with -g to get the list of 
                   valid strains
    -s, --strains  the list of strains to use, seperated by ':'
    -t, --tab      the probe file is tab delimited, the default for -f and -z
    -v, --verbose  show informational messages
    --vcf          the gzipsnp file is in vcf format.  if this is not used the
                   format is a format defined within the CGD, described below.
    -z, --zip      a zip containing the probe file. This also assumes there is a
                   file in the zip named probes.tsv, and the file is --tab
    
  Copyright (c) 2019 The Churchill Lab - The Jackson Laboratory
  
  This is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This software is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this software.  If not, see <http://www.gnu.org/licenses/>.

"""

from io import StringIO

import logging
import re
import sys

from . import exceptions
from . import desnp_utils



LOG = logging.getLogger('DeSNP')

#
#
# Logging
#
#


class DeSNPFormatter(logging.Formatter):

    dbg_fmt = "[DeSNP debug] %(msg)s"
    info_fmt = "[DeSNP] %(msg)s"
    err_fmt = "[DeSNP] %(msg)s"
    # err_fmt = "DEBUG: %(asctime)s, %(module)s: %(lineno)d: %(msg)s"
    # dbg_fmt = "DEBUG: %(asctime)s, %(module)s: %(lineno)d: %(msg)s"

    def __init__(self, fmt="[DeSNP] %(msg)s"):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._fmt = DeSNPFormatter.dbg_fmt

        elif record.levelno == logging.INFO:
            self._fmt = DeSNPFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._fmt = DeSNPFormatter.err_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._fmt = format_orig

        return result


def get_logger():
    """
    Get the logger
    :return: logger object
    """
    return LOG


def configure_logging(level):
    """

    :param level: 0=ERROR, 1=INFO, 2+=DEBUG
    """
    global LOG
    LOG = logging.getLogger('DeSNP')

    handler = logging.StreamHandler(sys.stderr)

    mfmt = DeSNPFormatter()
    handler.setFormatter(mfmt)

    if level == 0:
        LOG.setLevel(logging.INFO)
    else:
        LOG.setLevel(logging.DEBUG)

    LOG.addHandler(handler)


def exit(message='', parser=None):
    """
    Print message, help, and exit.

    :param message: The message to print
    :param parser: the argument parser
    """
    if parser:
        parser.error(str(message))
    else:
        if message:
            sys.stderr.write("[DeSNP] Error: {}\n".format(str(message)))

    sys.exit(1)




# The first column in the Sanger VCF where there are strains
VCF_STRAIN_START_COL = 9
VCF_SNP_POS = 1

# The first column in Dan's CGD SNP file where there are strains
# mm9 is 5, mm10 is 4... we need to do something about this.
CGD_STRAIN_START_COL = 5
#CGD_STRAIN_START_COL = 4


# Same problem here... mm10 different than mm9
CGD_SNP_POS = 2
#CGD_SNP_POS = 1
CGD_REF_POS = 3
#CGD_REF_POS = 2
CGD_ALT_POS = 4
#CGD_ALT_POS = 3


PROBE_FILE = "probes.tsv"
FILTERED_PROBE_FILE = "probes_filtered.tsv"
SNP_PROBE_FILE = "probes_snp.tsv"

# PROBE ID COLUMN IS REQUIRED.  We assume it is unique and named "id"
PROBE_ID_COL_NAME = "id"

#  The list of chromosomes supported for DeSNPing
CHRS = [str(x) for x in range(1, 20)]
CHRS.append('X')


def _deSNP(probe, probe_snp_list, strains, probe_writer, rej_writer, vcf=False, flag=False):
    """
    deSNP  Take a probes list of snps and sees if any are in the set of strains

    This function takes the probe,
    the list of snps found in the region of a probe,
    the list of strains of interest,
    the writer for outputing Probes,
    the writer for outputing probes rejected with a variance
    and a boolean whether or not the snps are in vcf format

    If a snp to the reference genome is found to exist in any of the strains
    of interest this entire probe is "DeSNPed", meaning dropped from the set
    of interest.  Probes that are DeSNPed are written to the rej_writer with an
    extra column containing the position and strain(s) that had the snp.  The
    probes with no SNPs between strains are written to the probe_writer.

    The "strain/SNP" column of the rejected file is formated:
       strain1;strain2;strain3;...;
     for the header and:
       0:1;1;;...;0:2
     where under each strain is the colon separated list of positions
     with a SNP for that strain, empty string where there are no SNPs
     for the strain.

    There is no return for this method.
    """
    # snpmap holds the list of snps found by strain
    snpmap = {}
    for strain in strains.keys():
        # Initialize with "", meaning no snps for this strain
        #
        snpmap[strain] = ""

    locations = []
    if probe.location:
        tokens = probe.location.split(";")
        for token in tokens:
            (chrom, prange) = token.split(":")
            (start, end) = prange.split("-")
            locations.append({"start": start, "end": end})
        if len(locations) > 1:
            # If negative strand, the locations are stored highest start
            # to lowest
            if probe.strand == "-":
                locations.reverse()
    else:
        locations.append({"start": probe.probe_start, "end": probe.probe_end})

    # assume there are no snps within the set of strains
    has_snp = False
    for snp in probe_snp_list:
        tokens = snp.split("\t")
        snp_position = -1
        for strain in sorted(strains.keys(), key=str.lower):
            # This is the position in the lookup row where strain is located
            position = strains[strain]
            # column1 in the tokens list is "SNP Position"
            if vcf:
                snp_position = tokens[VCF_SNP_POS]
            else:
                snp_position = tokens[CGD_SNP_POS]
            ref = ""
            alt = ""
            confidence = ""
            if vcf:
                #  Adjust the position based on start of strains in row
                position = position + VCF_STRAIN_START_COL
            else:
                # reference call
                ref = tokens[CGD_REF_POS]
                # alternate call
                alt = tokens[CGD_ALT_POS]
                #  Adjust the position based on start of strains in row
                position = ((position) * 2) + CGD_STRAIN_START_COL
                confidence = tokens[position + 1]

            value = tokens[position]
            #  if in the case of VCF
            if ((vcf and (value.startswith("1/1") or value.startswith("0/1") or
                value.startswith("1/0"))) or
                (value == alt and (confidence == "1" or confidence == "2"))):
                has_snp = True
                corrected_position = 0
                if (len(locations) > 1):

                    corrected_position = int(snp_position)
                    last_end = 0
                    for location in locations:
                        if snp_position >= location["start"]:
                            if last_end == 0:
                                corrected_position -= int(location["start"])
                                last_end = int(location["end"])
                            else:
                                diff = int(location["start"]) - last_end
                                corrected_position -= diff
                                last_end = int(location["end"])
                        else:
                            break
                else:
                    #  TODO:  This shouldn't work!!!
                    corrected_position = int(snp_position) - int(locations[0]["start"])
                # If no snp found yet, replace '0' with single letter abbrev
                if snpmap[strain] == '':
                    snpmap[strain] = str(corrected_position)
                    #  This is a debug version of above line that writes with strain name
                    #snpmap[strain] = strain + "-" + str(corrected_position)
                # If snp already found here, append single letter abbrev
                else:
                    snpmap[strain] = snpmap[strain] + ":" + str(corrected_position)
                    #  This is a debug version of above line that writes with strain name
                    #snpmap[strain] = snpmap[strain] + ":" + strain + "-" + str(corrected_position)
            #else:
            #    if snpmap[strain] == '':
            #        snpmap[strain] == '#'
            #
            #    else:
            #        snpmap[strain] = snpmap[strain] + ":" + '#'

    if has_snp:
        # The "strain/SNP" column of the reject file is formated:
        #   strain1;strain2;strain3;...;strainN
        # for the header and:
        #   0:1;1;;...;0:2
        # where under each strain is the colon separated list of position
        # with a SNP for that strain, empty string where there are no SNPs
        # for the strain.
        # Concatenate snp positions
        snpstring = ""
        first = True
        for strain in sorted(snpmap.keys(), key=str.lower):
            if first:
                snpstring = str(snpmap[strain])
                first = False
            else:
                snpstring = snpstring + ";" + str(snpmap[strain])
        probe_list = probe.asList()
        probe_list.append(snpstring)
        # write to rejected variance file
        rej_writer.writerow(probe_list)
    else:
        # write to probe file
        out_list = probe.asList()
        if flag:
            out_list.append("")
        probe_writer.writerow(out_list)

    return (1, 0) if has_snp else (0, 1)


def get_strains(snps_file, is_vcf=False):
    strains = []

    try:
        snps_file = desnp_utils.check_file(snps_file)

        if is_vcf:
            strains = get_samples_vcf(snps_file)
        else:
            strains = get_samples_cgd(snps_file)
    except Exception as e:
        raise e

    return strains


def get_SNP_list(probes, tabix_file):
    """
    get_SNP_list()  returns the list of snps found in probe region

    This function takes a probe, and a snp file.

    The function returns a list of snps, if any are found in the region
    covered by this probe.  If no snps are found "None" is returned.
    """
    regions = []
    for probe in probes:
        if probe.probe_start > -1:
            tmp_regions = tabix_file.fetch(reference=probe.chromosome,
                                           start=probe.probe_start,
                                           end=probe.probe_end)
            regions += tmp_regions
    return regions


def get_samples_vcf(snp_file_name):
    """
    get_samples_vcf()  Gets the list of valid sample/strains from a VCF file

    Takes the name of the gzipped vcf file as a parameter.  Assumes that in
    the VCF format, that the first line with only one '#' is the header line.
    Also assumes that the 9th column is the first one with a sample name.

    Returns a list of  strain names
    """
    with gzip.open(snp_file_name, 'rt') as gz:
        line = gz.readline()
        strains = []

        while line is not None:
            if line.startswith("##"):
                line = gz.readline()
                continue
            elif line.startswith("#"):
                tokens = line.split("\t")
                for col in range(VCF_STRAIN_START_COL, len(tokens)):
                    strains.append(tokens[col].strip())
                break
            else:
                break

        LOG.debug("Found {} strains ... {}".format(len(strains), str(strains)))
    return strains


def get_samples_cgd(snp_file_name):
    """
    get_samples_cgd()  Gets list of valid sample/strains from CGD strain file

    Takes the name of the gzipped cgd strain file as a parameter.  Assumes that in
    the CGD Strain format, the first line is a header.
    Also assumes that the 5th column is the first one with a sample name, and
    that every other column after that is a strain and that the alternating
    columns are a confidence score.

    Returns a list of  strain names
    """
    strains = []

    with gzip.open(snp_file_name, 'rt') as gz:
        line = gz.readline()

        # The first line is the header, split it...

        tokens = line.split("\t")

        # Once we start parsing strains every other column is a strain,
        # The others are a confidence.  For now we are not keeping the conf val

        strain = True
        for col in range(CGD_STRAIN_START_COL, len(tokens)):
            if strain:
                strains.append(tokens[col].strip())
                strain = False
            else:
                strain = True

        LOG.debug("Found {} strains ... {}".format(len(strains), str(strains)))

    return strains


def validate_strains(strains_string, snp_file_name, vcf=False):
    """
    validate_strains()  validates strains of interest against snp file.

    Takes a string of strain names separated by ":", a gzipped snp file, and
    a boolean flag as to whether or not that file is VCF format.

    Checks to see if the strains requested are avaiable in the snp file.

    If they are, a dict of strains is returned with the key=strain name and
    value=column in snp file.

    If any are not, an error is thrown with appropriate message and program
    exits.
    """
    if vcf:
        valid_samples = get_samples_vcf(snp_file_name)
    else:
        valid_samples = get_samples_cgd(snp_file_name)

    strains = strains_string.split(":")
    strain_dict = {}

    for strain in strains:
        try:
            strain_dict[strain] = valid_samples.index(strain)
        except ValueError:
            msg = "Strain '{}' not in valid set.".format(strain)
            LOG.error(msg)
            LOG.error("Valid strains available:")
            LOG.error("------------------------")

            for strain in valid_samples:
                LOG.error("    {}".format(strain))

            raise exceptions.DeSNPError(msg)

    return strain_dict


def perform_desnp(probes_file, snps_file, strains_string, output_file,
                  delimiter='\t', is_vcf=False, id_col="", flag=False):
    """

    :param probes_file:
    :param snps_file:
    :param strains_string:
    :param output_file:
    :param delimiter:
    :param is_vcf:
    :param id_col:
    :param flag:
    :return:
    """
    try:
        probes_file = desnp_utils.check_file(probes_file)
        snps_file = desnp_utils.check_file(snps_file)
        output_file = desnp_utils.check_file(output_file, mode="w")

        strains = validate_strains(strains_string, snps_file, is_vcf)

        with desnp_utils.open_resource(probes_file, 'rt') as probes_fd, \
             desnp_utils.open_resource(output_file, 'w') as out_fd:

            reader = csv.reader(probes_fd, delimiter=delimiter)
            writer = csv.writer(out_fd, delimiter=delimiter)

            # Set up writer for probes with snps (rejected)
            # This file will contain our probes that were rejected because of
            # variances/SNPs in between our selected strains
            rej_fd = None

            if not flag:
                rej_file_name = SNP_PROBE_FILE
                rej_fd = open(rej_file_name, 'w')
                rej_writer = csv.writer(rej_fd, delimiter=delimiter)
            else:
                rej_writer = writer

            snp_reader = pysam.Tabixfile(snps_file, 'r')

            first = True
            probe_counter = 0
            written_probes = 0
            written_snps = 0
            input_header = None
            headers_written = False

            DEBUG_TIME_TEST = {"findSnp": None, "deSnp": None}

            # process the probe file
            for line in reader:
                # the first line should be a header, write it back out...
                if first:
                    first = False
                    input_header = line
                    continue

                tmp_probes = parse_probe(line, input_header, PROBE_ID_COL_NAME)
                probe_counter += len(tmp_probes)

                # if we haven't written out the header line for our output files
                # do it now.  Must be done after we have our first probe.
                if not headers_written:
                    # write header to filtered file
                    head_list = tmp_probes[0].headList()

                    if not flag:
                        writer.writerow(head_list)

                    rej_head = head_list

                    # The "strain/SNP" column of the rejected file is formatted:
                    #   strain1;strain2;strain3;...;strainN
                    # for the header and:
                    #   0:1;1;;...;0:2
                    # where under each strain is the colon separated list of
                    # position with a SNP for that strain, empty string where
                    # there are no SNPs for the strain.
                    rej_strains = ";".join(sorted(strains.keys(), key=str.lower))
                    rej_head.append(rej_strains)
                    rej_writer.writerow(rej_head)
                    headers_written = True

                # Checking to see if this is a valid Chromosome...ValueError would be a no.
                try:
                    #  All probes in tmp_probes are the same probe, different locations, so the
                    #  chromosome is the same
                    CHRS.index(tmp_probes[0].chromosome)
                except ValueError:
                    if tmp_probes[0].chromosome is None:
                        continue
                    else:
                        #  Chromosome not supported, keep it in probe set...
                        # logging.warning("Chr " + str(tmp_probes[0].chromosome) + " not supported.  " +
                        #                "Keeping probe... " + str(tmp_probes[0].id))
                        out_list = tmp_probes[0].asList()
                        if flag:
                            out_list.append("")
                        writer.writerow(out_list)
                        written_probes += 1
                        continue
                else:
                    a = datetime.now()
                    probe_snp_list = get_SNP_list(tmp_probes, snp_reader)
                    b = datetime.now()
                    c = b - a
                    if DEBUG_TIME_TEST["findSnp"] == None:
                        DEBUG_TIME_TEST["findSnp"] = c
                    else:
                        DEBUG_TIME_TEST["findSnp"] = DEBUG_TIME_TEST["findSnp"] + c

                    a = datetime.now()
                    ret = _deSNP(tmp_probes[0], probe_snp_list, strains, writer, rej_writer, is_vcf, flag)
                    written_snps += ret[0]
                    written_probes += ret[1]
                    b = datetime.now()
                    c = b - a

                    if DEBUG_TIME_TEST["deSnp"] == None:
                        DEBUG_TIME_TEST["deSnp"] = c
                    else:
                        DEBUG_TIME_TEST["deSnp"] = DEBUG_TIME_TEST["deSnp"] + c

            LOG.debug("Took " + str(DEBUG_TIME_TEST["findSnp"]) + " to find SNPS")
            LOG.debug("Took " + str(DEBUG_TIME_TEST["deSnp"]) + " to do DeSNPing")
            LOG.info("finished processing at: " + time.strftime("%H:%M:%S") +
                         "\n")
            LOG.info("Total probes = " + str(probe_counter) + \
                         " probes without snps between strains = " + str(written_probes) \
                         + " probes with snps between strains = " + str(written_snps))

            LOG.info("DeSNP Completed Successfully")
            sys.exit(0)
    except Exception as detail:
        #import traceback
        #traceback.print_exc(file=sys.stdout)
        #print("-" * 60)
        LOG.error("DeSNP did not run to completion " + str(detail))
        sys.exit(1)

