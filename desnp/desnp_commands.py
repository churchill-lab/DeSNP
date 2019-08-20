# -*- coding: utf-8 -*-

import argparse
import sys

from . import desnp
from . import exceptions
from . import summarize


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
    parser.add_argument("-o", "--out", dest="output", metavar="Output_File",
                        default="probes.out")
    parser.add_argument("--vcf", dest="vcf", action='store_true')

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count",
                        default=0)

    args = parser.parse_args(raw_args)

    desnp.configure_logging(args.debug)

    if args.help:
        desnp.exit("", parser)

    print(args)

    try:
        if not args.snps:
            desnp.exit("No snps file was specified.", parser)

        if not args.strains:
            desnp.exit("No strains were specified.", parser)

        if not args.probes:
            desnp.exit("No probes file was specified.", parser)

        delimiter = ',' if args.comma else '\t'

        desnp.perform_desnp(args.probes, args.snps, args.strains, args.output,
                            delimiter, args.vcf, args.idcol, args.flag)

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

    Usage: DeSNP strains <snpsfile.gz>

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
    parser.add_argument("-d", "--debug", dest="debug", action="count",
                        default=0)

    args = parser.parse_args(raw_args)

    desnp.configure_logging(args.debug)

    if args.help:
        desnp.exit("", parser)

    try:
        if not args.snpsfile:
            desnp.exit("No snps file was specified.", parser)

        strains = desnp.get_strains(args.snpsfile, args.vcf)

        if strains == 1:
            LOG.error("Please provide a SNP File.")
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
    Usage: DeSNP summarize -f <data.tsv> -p <probes.tsv> -s <samples.tsv>

    Required Parameters:

        -f, --data <data file>           The matix of intensity data

        -p, --probes <probes file>       Input probes file, tab delimited.

        -s, --samples <samples file>     The design file containing the samples


    Optional Parameters:

        -c, --samplecol                  The column in the design file
                                         containing the unique sample
                                         identifier/name

        -e, --extra                      Generate an additional json file
                                         containing extra median polish results.
                                         This option only works with gene
                                         grouping, otherwise median polish not
                                         run

        -g, --group                      How to group probe sets, options are:
                                         probe, gene (default)

        -i, --idcol                      The name of the unique probe id column.
                                         If not provided assumes 'id'

        -o, --out                        The name of the output file the results
                                         will go to


    Help Parameters:

        -h, --help                       Print the help and exit

        -d, --debug                      Turn debugging mode on (list multiple
                                         times for more messages)

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
    parser.add_argument("-f", "--data", dest="data", metavar="Data_File")
    parser.add_argument("-p", "--probes", dest="probes", metavar="Probes_File")
    parser.add_argument("-s", "--samples", dest="samples",
                        metavar="Samples_File")

    # optional
    parser.add_argument("-c", "--samplecol", dest="samplecol",
                        metavar="Sample_Column")
    parser.add_argument("-e", "--extra", dest="extra", action='store_true')
    parser.add_argument("-g", "--group", dest="group",
                        choices=['probe', 'gene'], default="gene")
    parser.add_argument("-i", "--idcol", dest="idcol", metavar="ID_Column")
    parser.add_argument("-o", "--out", dest="output", metavar="Output_File",
                        default="results.out")

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count",
                        default=0)

    args = parser.parse_args(raw_args)

    desnp.configure_logging(args.debug)

    if args.help:
        desnp.exit("", parser)

    try:
        if not args.data:
            desnp.exit("No data file was specified.", parser)

        if not args.probes:
            desnp.exit("No probes file was specified.", parser)

        if not args.samples:
            desnp.exit("No sample file was specified.", parser)

        summary = summarize.Summary(args.data, args.probes, args.samples,
                                    args.output, args.group, args.extra,
                                    args.idcol, args.samplecol)

        summary.process()
        LOG.info("Summarization Completed Successfully")

    except KeyboardInterrupt as ki:
        LOG.warn(ki)
    except exceptions.DeSNPError as ve:
        desnp.exit(ve, parser)
    except Exception as detail:
        import traceback
        print("-" * 60)
        traceback.print_exc(file=sys.stdout)
        LOG.info("Summarization did not run to completion " + str(detail))
        sys.exit(1)
