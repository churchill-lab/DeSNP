# -*- coding: utf-8 -*-

import argparse
import sys

from . import desnp
from . import exceptions

LOG = desnp.get_logger()


def command_process(raw_args, prog=None):
    """
    DeSNP

    Usage: DeSNP process -f <probes.txt> -g <snps.gz> -s <strains>

    Required Parameters:

        -p, --probes <probes file>       Input probes file

        -s, --snps <snps file>           Gzipped snp file.  This requires an
                                         associated .tbi tabix index file to
                                         be present in the same location

        -t, --strains <strains>          List of strains, separated by ':'


    Optional Parameters:

        -c, --comma                      The probe file is comma delimited,
                                         defaults to TAB delimited.

        --flag                           Default behavior of DeSNP is to filter
                                         out probes with snps and write them to
                                         a separate file.  The --flag option
                                         keeps them in one probe file with an
                                         additional column added to identify
                                         location and strain of snps

        -i, --idcol                      The name of the unique probe id column.
                                         If not provided assumes 'id'

        -o, --out                        The name of the output file the results
                                         will go to.

        --vcf                            The gzipsnp file is in vcf format.  If
                                         this is not used the format is a format
                                         defined within the CGD, described below

    Help Parameters:

        -h, --help                       Print the help and exit

        -d, --debug                      Turn debugging mode on (list multiple
                                         times for more messages)


    Note:
            CGD SNP File format: Tab delimited file containing the following
            columns - SNPID, CHROM, POS, REF, ALT, Strain1 Allele, Strain1
            confidence, ... StrainN Allele, StrainN conf
    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        if message:
            sys.stderr.write(message)
        else:
            sys.stderr.write(command_process.__doc__)
        sys.stderr.write('\n')
        sys.exit(1)

    parser.error = print_message

    # required
    parser.add_argument("-p", "--probes", dest="probes", metavar="Probes_File")
    parser.add_argument("-s", "--snps", dest="snps", metavar="VCI_File")
    parser.add_argument("-t", "--strains", dest="strains", metavar="Strains")

    # optional
    parser.add_argument("-c", "--comma", dest="comma", action='store_true')
    parser.add_argument("--flag", dest="flag", action='store_true')
    parser.add_argument("-i", "--idcol", dest="idcol", metavar="ID_Column")
    parser.add_argument("-o", "--out", dest="output", metavar="Output_File", default="probes.out")
    parser.add_argument("--vcf", dest="vcf", action='store_true')

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    desnp.configure_logging(args.debug)

    if args.help:
        desnp.exit("", parser)

    try:
        if not args.snps:
            desnp.exit("No snps file was specified.", parser)

        if not args.strains:
            desnp.exit("No strains were specified.", parser)

        if not args.probes:
            desnp.exit("No probes file was specified.", parser)

        delim = ',' if args.comma else '\t'

        desnp.perform_desnp(args.probes, args.snps, args.strains, args.output,
                            delim, args.vcf, args.idcol, args.flag)

    except KeyboardInterrupt as ki:
        LOG.warn(ki)
    except exceptions.DeSNPError as ve:
        desnp.exit(ve, parser)
    except Exception as detail:
        LOG.error("DeSNP did not run to completion {}".format(str(detail)))
        sys.exit(1)


def command_strains(raw_args, prog=None):
    """
    DeSNP

    Usage: DeSNP strains <snps.gz>

    Required Parameters:

        <snps file>                      Gzipped snp file.  This requires an
                                         associated .tbi tabix index file to
                                         be present in the same location

    Optional Parameters:

        --vcf                            The gzipsnp file is in vcf format.  If
                                         this is not used the format is a format
                                         defined within the CGD, described below

    Help Parameters:

        -h, --help                       Print the help and exit

        -d, --debug                      Turn debugging mode on (list multiple
                                         times for more messages)


    Note:
            CGD SNP File format: Tab delimited file containing the following
            columns - SNPID, CHROM, POS, REF, ALT, Strain1 Allele, Strain1
            confidence, ... StrainN Allele, StrainN conf
    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        if message:
            sys.stderr.write(message)
        else:
            sys.stderr.write(command_strains.__doc__)
        sys.stderr.write('\n')
        sys.exit(1)

    parser.error = print_message

    # required
    parser.add_argument("snpsfile", metavar="snpsfile.gz")

    # optional
    parser.add_argument("--vcf", dest="vcf", action='store_true')

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    desnp.configure_logging(args.debug)

    if args.help:
        desnp.exit("", parser)

    try:
        if not args.snpsfile:
            desnp.exit("No snps file was specified.", parser)

        strains = desnp.get_strains(args.snpsfile, args.vcf)

        if strains == 1:
            print("Please provide a SNP File.")
            sys.exit(1)
        else:
            LOG.info("Valid strains available:")
            LOG.info("------------------------")

            for strain in strains:
                LOG.info("    {}".format(strain))

            sys.exit(0)
    except KeyboardInterrupt as ki:
        LOG.warn(ki)
    except exceptions.DeSNPError as ve:
        desnp.exit(ve, parser)
    except Exception as detail:
        LOG.error("DeSNP did not run to completion {}".format(str(detail)))
        sys.exit(1)


def command_summarize(raw_args, prog=None):
    """
    Create VCI file from VCF file(s)

    Usage: g2gtools vcf2vci [-i <indel VCF file>]* -s <strain> -o <output file> [options]

    Required Parameters:
        -i, --vcf <vcf_file>             VCF file name
        -f, --fasta <Fasta File>         Fasta file matching VCF information
        -s, --strain <Strain>            Name of strain (column in VCF file)
        -o, --output <Output file>       VCI file name to create

    Optional Parameters:
        -p, --num-processes <number>     The number of processes to use, defaults to the number of cores
        --diploid                        Create diploid VCI
        --keep                           Keep track of VCF lines that could not be converted to VCI file
        --pass                           Use only VCF lines that have a PASS for the filter value
        --quality                        Filter on quality, FI=PASS
        --no-bgzip                       DO NOT compress and index output

    Help Parameters:
        -h, --help                       Print the help and exit
        -d, --debug                      Turn debugging mode on (list multiple times for more messages)

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        if message:
            sys.stderr.write(message)
        else:
            sys.stderr.write(command_summarize.__doc__)
        sys.stderr.write('\n')
        sys.exit(1)

    parser.error = print_message

    # required
    parser.add_argument("-i", "--vcf", dest="vcf_files", metavar="vcf_file", action='append')
    parser.add_argument("-f", "--fasta", dest="fasta_file", metavar="fasta_file")
    parser.add_argument("-s", "--strain", dest="strain", metavar="strain")
    parser.add_argument("-o", "--output", dest="output", metavar="VCI_File")

    # optional
    parser.add_argument("-p", "--num-processes", type=int, dest="numprocesses", metavar="number_of_processes")
    parser.add_argument("--diploid", dest="diploid", action='store_true')
    parser.add_argument("--keep", dest="keep", action='store_true')
    parser.add_argument("--pass", dest="passed", action='store_true')
    parser.add_argument("--quality", dest="quality", action='store_true')
    parser.add_argument("--no-bgzip", dest="nobgzip", action='store_false', default=True)

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    desnp.configure_logging(args.debug)


