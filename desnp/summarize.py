# -*- coding: utf-8 -*-

import csv
import json
import operator
import time
from datetime import datetime

import numpy as np

from . import desnp
from . import desnp_utils
from . import exceptions
from . import medpolish as mp
from . import probe

LOG = desnp.get_logger()


def quantnorm(x):
    """
    quantnorm() is a method for doing quantile normalization.
    It takes a numpy matrix and does quantile normalization of the matrix.
    Returns the normalized matrix
    """
    rows, cols = x.shape
    sortIndexes = np.argsort(x, axis=0)

    for row in range(rows):
        currIndexes = sortIndexes[row, :]
        x[currIndexes, range(cols)] = np.median(x[currIndexes, range(cols)])


class Summary(object):

    def __init__(self, data_file, probe_file, sample_file,
                 out_file=None, group='gene', extra=False,
                 idcol=None, sample_col=None):

        #  This list is used to keep the order we read our probes in.
        self.probe_ids = []

        #  This indicates whether we are grouping by probe or gene
        self.group = None

        self.group = group
        self.extra = extra
        self.delim = '\t'
        self.probe_id_col_name = "id"
        self.probe_set_col_name = "ProbeSet ID"
        self.sample_col_name = "sampleid"
        self.sample_col_name_alt = "Sample"

        self.data_file = desnp_utils.check_file(data_file)
        self.probe_file = desnp_utils.check_file(probe_file)
        self.sample_file = desnp_utils.check_file(sample_file)

        if out_file is not None:
            self.out_file = desnp_utils.check_file(out_file, 'w')
        else:
            self.out_file = desnp_utils.check_file('statistics.tsv', 'w')

        if idcol is not None:
            self.probe_id_col_name = idcol

        if sample_col is not None:
            self.sample_col_name = sample_col

        LOG.debug("Using data file provided: {}".format(data_file))
        LOG.debug("Using probe file provided: {}".format(probe_file))
        LOG.debug("Using sample file provided: {}".format(sample_file))

    def process(self):

        LOG.debug("Started processing at: {}".format(time.strftime("%H:%M:%S")))

        with desnp_utils.open_resource(self.data_file, 'rt') as data_fd, \
             desnp_utils.open_resource(self.probe_file, 'rt') as probe_fd, \
             desnp_utils.open_resource(self.sample_file, 'rt') as sample_fd, \
             desnp_utils.open_resource(self.out_file, 'w') as writer_fd:

            writer = csv.writer(writer_fd, delimiter=self.delim)

            _time = datetime.now()

            # get our set of filtered probes
            probes = self.get_probes(probe_fd)

            # get our set of sample names
            samples = self.get_sample_names(sample_fd)

            # Add the intensity data to the probe objects
            self.add_probe_data(probes, data_fd)

            LOG.debug("Loading data took: {}".format(datetime.now() - _time))

            # diagnostic count only
            group_count = 0

            # if there are no probes left after filtering for gene start and end
            # skip this step
            if len(probes) > 0:
                # Generate a single intensity matrix
                intensity_matrix = self.get_intensity_matrix(probes)

                _time = datetime.now()

                # log2 Transform matrix
                log2_matrix = np.log2(intensity_matrix)

                LOG.debug("log2 took: {}".format(datetime.now() - _time))

                _time = datetime.now()

                # do quantile norm of matrix (method updates log2_matrix by ref)
                quantnorm(log2_matrix)

                LOG.debug("quantnorm took: {}".format(datetime.now() - _time))

                i = 0

                _time = datetime.now()

                for probe_id in self.probe_ids:
                    intensity_row = log2_matrix[i]
                    inten_array = intensity_row.tolist()
                    probe = probes[probe_id]
                    probe.setIntensities(inten_array)
                    i += 1

                LOG.debug("resetting intensities took: {}".format(
                    datetime.now() - _time))

                # create our groupings and write the the statistics file
                if self.group != 'probe':
                    _time = datetime.now()
                    groupings = self.group_probesets_by_gene(probes, samples)

                    LOG.debug("grouping Probesets took: {}".format(
                        datetime.now() - _time))

                    group_count = len(groupings)

                    _time = datetime.now()
                    sorted(groupings, key=lambda grouping: (grouping.chromosome,
                                                            grouping.start_pos))

                    LOG.debug("sorting groupings took: {}".format(
                        datetime.now() - _time))

                    first = True

                    _time = datetime.now()
                    medpol_groupings = []

                    for grouping in groupings:
                        if first:
                            writer.writerow(grouping.headList())
                            first = False
                        writer.writerow(grouping.asList())

                        if self.extra:
                            medpol_groupings.append(grouping.med_polish_results)

                    if self.extra:
                        md = open("median_polish_by_group.json", 'w')
                        md.write(json.dumps(medpol_groupings))

                    LOG.debug("writing groupings took: {}".format(
                        datetime.now() - _time))
                else:
                    first = True
                    group_count = len(self.probe_ids)

                    for probe_id in self.probe_ids:
                        if first:
                            writer.writerow(probe.headList())
                            first = False

                        writer.writerow(probe.asList())
            else:
                LOG.info("No probes left to summarize.")
                LOG.info("Possibly there were no Gene Start and End values in")
                LOG.info("'{}'".format(self.probe_file))

            # force any data written to file
            writer_fd.flush()
            writer_fd.close()

        LOG.info("finished processing at: {}".format(time.strftime("%H:%M:%S")))
        LOG.info("Total groupings = {} from {} probes".format(
            str(group_count), str(len(self.probe_ids))))

    def get_probes(self, probe_fd):
        """
        Takes a file descriptor for a probe file.  Assumes file is delimited by
        tabs. Makes the assumption that the probe set id column is "ProbeSet ID"
        and that the probe id column is "id".  Does not assume column order,
        but uses these two fields to extract information.  Returns a dictionary
        of Probe objects.
        """
        probes = {}
        reader = csv.reader(probe_fd, delimiter="\t")
        first = True
        header = []
        id_col = -1
        psi_col = -1
        dropped_no_start_end = 0

        for line in reader:
            # skip header
            if first:
                header = line
                LOG.debug("Header: {}".format(header))

                if self.probe_set_col_name in header:
                    psi_col = header.index(self.probe_set_col_name)

                if self.probe_id_col_name in header:
                    id_col = header.index(self.probe_id_col_name)
                    LOG.debug("ID COL IS: {}".format(id_col))
                else:
                    msg = "ID column name {} does'nt exist in probe input file!"
                    msg = msg.format(self.probe_id_col_name)
                    LOG.error(msg)
                    raise exceptions.DeSNPSummarizationError(msg)

                first = False
                continue

            _probe = None

            #  a probe can appear many times in the file due to different
            #  probe set ids, all other columns are the same, so we just
            #  add the probesetid

            probe_id = line[id_col]

            if probe_id in probes and psi_col > -1:
                _probe = probes[probe_id]
                _probe.addProbeSetId(line[psi_col])
            else:
                # parse_prob can result in multiple instances of the same probe
                tmp_probes = \
                    probe.parse_probe(line, header,
                                      probe_id_col_name=self.probe_id_col_name)

                # for summarization this is irrelevant, just keep the first
                _probe = tmp_probes[0]

                # If grouping by gene and no gene start and end position, than
                # drop this probe
                if self.group == "gene" and \
                        (not _probe.start_pos or not _probe.end_pos):
                    dropped_no_start_end += 1
                    continue

                probes[_probe.id] = _probe
                self.probe_ids.append(_probe.id)

        LOG.info("Loaded {} probes.".format(str(len(self.probe_ids))))
        LOG.info("{} probes dropped.".format(str(dropped_no_start_end)))

        return probes

    def get_sample_names(self, sample_fd):
        """
        Takes the sample file descriptor and returns a list of sample names.
        It is assumed the sample names are in the column "sampleid" or "Sample"
        """
        reader = csv.reader(sample_fd, delimiter="\t")
        samples = []
        first = True
        sample_col = 0

        for line in reader:
            if first:
                try:
                    sample_col = line.index(self.sample_col_name)
                except ValueError:
                    try:
                        sample_col = line.index(self.sample_col_name_alt)
                    except ValueError:
                        # Fail, there required column name is missing
                        msg = """The sample file, must contain a column '{}' or
                                '{}'. This column should contain the names for 
                                the sample column headers in the summarized 
                                output file.""".format(self.sample_col_name,
                                                       self.sample_col_name_alt)
                        LOG.error(msg)
                        raise exceptions.DeSNPSummarizationError(msg)
                first = False
            else:
                samples.append(line[sample_col])

        LOG.debug("Sample names: ")
        LOG.debug(str(samples))

        return samples

    def add_probe_data(self, probes, data_fd):
        """
        Takes the map of probe objects and a file descriptor for the matrix of
        intensity values and adds these intensity values to the appropriate
        probe.
        """
        reader = csv.reader(data_fd, delimiter="\t")
        first = True
        updated = 0
        lines_read = 0

        for line in reader:
            lines_read += 1

            # skip header
            if first:
                first = False
                continue

            # skip blank lines
            if len(line) == 0:
                continue

            if len(line) == 1:
                msg = "Only 1 column on line {}.  May be using wrong delimiter."
                msg = msg.format(lines_read)
                LOG.error(msg)
                raise exceptions.DeSNPSummarizationError(msg)

            probe_id = line[0]

            try:
                probes[probe_id].setIntensities(line[1:])
                updated += 1
            except:
                # If the probe_id is not in the dictionary, we skip the line
                continue

        LOG.debug("Updated {} probes with intensity data".format(str(updated)))

    def get_intensity_matrix(self, probes):
        """
        Will take the dictionary of probes and will return a matrix of
        intensity values probes x samples
        """
        data = []
        row_length = 0
        row_num = 0
        success = True

        # Iterate through probes and make sure each row has the same number of
        # intensities. If not, warn the user and exit, we cannot proceed with
        # summarization.
        for probe_id in self.probe_ids:
            row = probes[probe_id].intensities

            if row_length == 0:
                row_length = len(row)
            elif row_length != len(row):
                probe = probes[probe_id]
                msg = "{} intensities found, {} expected for probe: {}, {}, {}"
                msg = msg.format(str(len(row)), str(row_length), probe.probe_id,
                                 probe.probeset_id, probe.sequence)
                LOG.error(msg)
                success = False

            data.append(probes[probe_id].intensities)
            row_num = row_num + 1

        if not success:
            msg = "Cannot proceed with summarization.  Exiting..."
            LOG.error(msg)
            raise exceptions.DeSNPSummarizationError(msg)

        x = np.array(data, dtype=np.float)
        y = np.array(x, dtype=np.int32)

        return y

    def group_probesets_by_gene(self, probes, samples):
        """
        Takes the Dictionary of probes and the list of samples and groups them
        by gene.  It uses median polish to give only one value per sample per
        grouping. The new grouped matrix is returned.
        """
        groupings = None

        # Currently we are only grouping by Gene.  If we add another level of
        # Grouping later we should either break out in another function, or
        # add an additional parameter here.

        groupings_dict = {}

        # Divide the dataset into groups by MGI ID
        no_gene_id = 0
        for probe_id in self.probe_ids:
            _probe = probes[probe_id]
            gene_id = _probe.gene_id

            if gene_id is None or gene_id == '':
                no_gene_id += 1
                continue

            if gene_id in groupings_dict:
                grouping = groupings_dict[gene_id]
                grouping.addProbe(_probe)
            else:
                group_values = [_probe.gene_id, _probe.symbol, _probe.name,
                                _probe.chromosome, _probe.start_pos,
                                _probe.end_pos, _probe.strand]
                group_header = ['Gene ID', 'Gene Symbol', 'Gene Name', 'Chr',
                                'Start', 'End', 'Strand']
                grouping = probe.ProbeSet(group_values, group_header)
                grouping.setSampleNames(samples)
                grouping.addProbe(_probe)

                groupings_dict[_probe.gene_id] = grouping

        LOG.debug("{} entries without a gene id".format(str(no_gene_id)))

        groupings = groupings_dict.values()

        LOG.debug("Have {} probesets...".format(str(len(groupings))))

        # for each group, take the set of probe intensity values and
        #    run them through medpolish, then add the "col" results as
        #    the intensity values of the group
        groups_processed = 0

        for grouping in groupings:
            matrix = grouping.getProbeNPMatrix()
            medp_result = mp.adjusted_medpolish(matrix)
            grouping.setIntensities(medp_result.col)
            grouping.setMedPolishResults(medp_result)
            groups_processed += 1

            if operator.mod(groups_processed, 1000) == 0:
                LOG.debug("Calculated median polish on {} groups at {}".format(
                    str(groups_processed), time.strftime("%H:%M:%S")))

        return groupings
