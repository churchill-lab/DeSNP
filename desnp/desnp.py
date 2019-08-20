# -*- coding: utf-8 -*-

import logging
import csv
import gzip
import sys
import time
from builtins import range

import pysam
from datetime import datetime
from desnp.probe import parse_probe


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
# CGD_STRAIN_START_COL = 4


# Same problem here... mm10 different than mm9
CGD_SNP_POS = 2
# CGD_SNP_POS = 1
CGD_REF_POS = 3
# CGD_REF_POS = 2
CGD_ALT_POS = 4
# CGD_ALT_POS = 3


PROBE_FILE = "probes.tsv"
FILTERED_PROBE_FILE = "probes_filtered.tsv"
SNP_PROBE_FILE = "probes_snp.tsv"

# PROBE ID COLUMN IS REQUIRED.  We assume it is unique and named "id"
PROBE_ID_COL_NAME = "id"

#  The list of chromosomes supported for DeSNPing
CHRS = [str(x) for x in range(1, 20)]
CHRS.append('X')


def _deSNP(probe, probe_snp_list, strains, probe_writer, rej_writer, vcf=False,
           flag=False):
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

                # adjust the position based on start of strains in row
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
                    corrected_position = \
                        int(snp_position) - int(locations[0]["start"])

                if snpmap[strain] == '':
                    # if no snp found yet, replace '0' with single letter abbrev
                    snpmap[strain] = str(corrected_position)
                else:
                    # if snp already found here, append single letter abbrev
                    snpmap[strain] = '{}:{}'.format(snpmap[strain],
                                                    str(corrected_position))

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

    Takes the name of the gzipped cgd strain file as a parameter.  Assumes that
    in the CGD Strain format, the first line is a header.
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
                  delimiter='\t', is_vcf=False, id_col=PROBE_ID_COL_NAME,
                  flag=False):
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

        if not id_col:
            id_col = PROBE_ID_COL_NAME

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

                tmp_probes = parse_probe(line, input_header, id_col)
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
                    rej_strains = ";".join(sorted(strains.keys(),
                                                  key=str.lower))
                    rej_head.append(rej_strains)
                    rej_writer.writerow(rej_head)
                    headers_written = True

                # checking to see if this is a valid Chromosome...
                # ValueError would be a no.
                try:
                    # all probes in tmp_probes are the same probe, different
                    # locations, so the chromosome is the same
                    CHRS.index(tmp_probes[0].chromosome)
                except ValueError:
                    if tmp_probes[0].chromosome is None:
                        continue
                    else:
                        # Chromosome not supported, keep it in probe set...
                        out_list = tmp_probes[0].asList()

                        if flag:
                            out_list.append("")

                        writer.writerow(out_list)
                        written_probes += 1
                        continue
                else:
                    _time = datetime.now()
                    probe_snp_list = get_SNP_list(tmp_probes, snp_reader)

                    if DEBUG_TIME_TEST["findSnp"] is None:
                        DEBUG_TIME_TEST["findSnp"] = datetime.now() - _time
                    else:
                        DEBUG_TIME_TEST["findSnp"] = \
                            DEBUG_TIME_TEST["findSnp"] + \
                            (datetime.now() - _time)

                    _time = datetime.now()
                    ret = _deSNP(tmp_probes[0], probe_snp_list, strains,
                                 writer, rej_writer, is_vcf, flag)

                    written_snps += ret[0]
                    written_probes += ret[1]

                    if DEBUG_TIME_TEST["deSnp"] is None:
                        DEBUG_TIME_TEST["deSnp"] = datetime.now() - _time
                    else:
                        DEBUG_TIME_TEST["deSnp"] = \
                            DEBUG_TIME_TEST["deSnp"] + \
                            (datetime.now() - _time)

            LOG.debug("Took {} to find SNPS".format(
                DEBUG_TIME_TEST["findSnp"]))

            LOG.debug("Took {} to do DeSNPing".format(
                DEBUG_TIME_TEST["deSnp"]))

            LOG.info("Finished processing at: {}".format(
                time.strftime("%H:%M:%S")))

            LOG.info("Total probes: {}".format(probe_counter))

            LOG.info("Probes without snps between strains: {}".format(
                written_probes))

            LOG.info("Probes with snps between strains: {}".format(
                written_snps))

            LOG.info("DeSNP Completed Successfully")
            sys.exit(0)
    except Exception as detail:
        import traceback
        traceback.print_exc(file=sys.stdout)
        print("-" * 60)
        LOG.error("DeSNP did not run to completion " + str(detail))
        sys.exit(1)

