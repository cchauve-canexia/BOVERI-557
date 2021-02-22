"""
Module handling alignments of reads clusters on amplicons
"""
from collections import defaultdict
from itertools import groupby
import logging
from textwrap import dedent

from Bio import pairwise2

from bin.clusters_utils import UNAMBIGUOUS
from bin.fastq_utils import phred_to_str
from bin.parameters_utils import (
    ALG_CLUSTER_MIN_SIZE,
    ALG_GE,
    ALG_GO,
    ALG_M,
    ALG_MM,
    ALG_NB_MAX,
    LOG_AMPLICON,
    LOG_SAMPLE,
)

# Gap symbol in alignments
GAP = str('-')

logger = logging.getLogger(__name__)

# Auxiliary functions

# CIGAR symbols
CIGAR_MATCH = '='  # Match
CIGAR_MISMATCH = 'X'  # Mismatch
CIGAR_INS = 'I'  # Insertion
CIGAR_DEL = 'D'  # Deletion
CIGAR_HARD_CLIP = 'H' # Hard clipping


def coord_to_aligned(aligned_seq1, aligned_seq2, pos):
    """
    Computes the position in the unaligned version of aligned_seq2 of the
    aligned position correspoding to pos in unaligned version of aligned_seq1
    :param: aligned_seq1, aligned_seq2 (str): aligned sequences
    :assumption: sequences of same length
    :param: pos (int): coordinate of a base in unaligne dversion of aligned_seq1
    :return: int
    """
    unaligned_seq1_pos, aligned_pos = 0, 0
    while unaligned_seq1_pos <= pos:
        unaligned_seq1_pos += (1 if aligned_seq1[aligned_pos] != GAP else 0)
        aligned_pos += 1
    aligned_pos -= 1
    return aligned_pos - aligned_seq2[0:aligned_pos].count(GAP)


def get_cigar_symbol(ref_pos, read_pos):
    """
    Get the cigar status of two aligned bases
    :param: ref_pos, read_pos (str): aligned bases
    :return: str: CIGAR symbol
    """
    if ref_pos == read_pos:
        return CIGAR_MATCH
    elif ref_pos == GAP:
        return CIGAR_INS
    elif read_pos == GAP:
        return CIGAR_DEL
    return CIGAR_MISMATCH


def compute_alignment_expanded_cigar(ref_seq, read_seq):
    """
    Computes expanded CIGAR string describing the alignment (ref, read)
    where each alignment position is matched to the corresponding CIGAR symbol
    :param: ref_seq, read_seq (str): aligned sequences
    (ref: reference sequence, read: read sequence)
    :return: str: expanded CIGAR string
    """
    # CIGAR symbol at every positions
    alg_pos = range(len(read_seq))
    return [get_cigar_symbol(ref_seq[pos], read_seq[pos]) for pos in alg_pos]


def compute_alignment_cigar(ref_seq, read_seq):
    """
    Computes CIGAR string describing the alignment (ref, read)
    :param: ref_seq, read_seq (str): aligned sequences
    (ref: reference sequence, read: read sequence)
    :return: str: CIGAR string
    """
    cigar_expanded = compute_alignment_expanded_cigar(ref_seq, read_seq)
    # Grouping consecutive identical CIGAR sybols
    return ''.join([f"{len(list(y))}{x}" for x, y in groupby(cigar_expanded)])


def compute_clipped_alignment_cigar(ref_aligned_seq,
                                    read_aligned_seq,
                                    read_clipped_seq,
                                    fwd=True):
    """
    Compute the CIGAR string of an alignment between a sequenced read
    (not a merged read), including the hard-clipped part
    :param: read_aligned_seq, ref_aligned_seq (str, str): aligned sequences
    :assumption: both sequences have the same length
    :param: read_clipped_seq (str): hard-clipped part of read
    :param: fwd (bool): True if the read is a forward read
    :return (str): CIGAR string
    """
    clipped_len = len(read_clipped_seq)
    clipped_cigar = ('' if clipped_len == 0 else f"{clipped_len}{CIGAR_HARD_CLIP}")
    if fwd:
        return (
            f"{compute_alignment_cigar(ref_aligned_seq, read_aligned_seq)}"
            f"{clipped_cigar}"
        )
    else:
        return (
            f"{clipped_cigar}"
            f"{compute_alignment_cigar(ref_aligned_seq, read_aligned_seq)}"
        )

def compute_alignment_NM(ref_seq, read_seq):
    """
    Computes the edit unit distance of an alignment (ref, read)
    :param: ref_seq, read_seq (str): aligned sequences
    (ref: reference sequence, read: read sequence)
    :return: int: edit distance
    """
    cigar_expanded = compute_alignment_expanded_cigar(ref_seq, read_seq)
    return (
        cigar_expanded.count(CIGAR_MISMATCH) +
        cigar_expanded.count(CIGAR_DEL) +
        cigar_expanded.count(CIGAR_INS)
    )


def compute_alignment_MD(ref_seq, read_seq):
    """
    Computes the MD tag string of an alignment (ref, read)
    :param: ref_seq, read_seq (str): aligned sequences
    (ref: reference sequence, read: read sequence)
    :return: int: MD tag string
    """
    # CIGAR symbol at every positions
    alg_pos = range(len(read_seq))
    MD, match_counter, current_state = '', 0, 'start'
    for pos in alg_pos:
        if ref_seq[pos] == read_seq[pos]: # match
            match_counter += 1
            current_state = 'match'
        elif read_seq[pos] != GAP and ref_seq[pos] != GAP: # mismatch
            MD += f"{match_counter}{ref_seq[pos]}"
            match_counter, current_state = 0, 'mismatch'
        elif read_seq[pos] == GAP: # Deletion
            if current_state != 'deletion':
                MD += f"{match_counter}^"
            MD += ref_seq[pos]
            match_counter, current_state = 0, 'deletion'
    MD += (str(match_counter) if match_counter > 0 else '')
    return MD

def compute_alignment_score(amp_seq, read_seq, alg_parameters):
    """
    Compute the score of an alignment
    :param: amp_seq, read_seq (str): aligned sequences
    :param: alg_parameters (dict(str, float)): alignment parameters
    :return: float: alignment score
    """

    def gap_opening(seq, pos):
        """ True if pos opens a gap in seq """
        return (seq[pos] == GAP and seq[pos - 1] != GAP)

    score = 0.0
    for pos in range(len(amp_seq)):
        amp, read = amp_seq[pos], read_seq[pos]
        if amp == read:
            score += alg_parameters[ALG_M]
        elif gap_opening(amp_seq, pos) or gap_opening(read_seq, pos):
            score += alg_parameters[ALG_GO]
        elif amp == GAP or read == GAP:
            score += alg_parameters[ALG_GE]
        else:
            score += alg_parameters[ALG_MM]
    return score

def align_single_read(alignment, seq_to_align, fwd=True):
    """
    Given an alignment of a cluster consensus sequence and a single read from
    this cluster where hard clipped parts have been removed, align the read to
    the reference sequence consistently to the initial cluster alignment
    :param: alignment (Alignment)
    :param: seq_to_align (str): alignable part of the read
    :param: fwd (bool): True if the read is a forward read
    :assumption: the aligned cluster sequence does not start or end with a gap
    :return: str, str: aligned read, realigned reference
    """
    ref_aligned_seq = alignment.get_amplicon()
    cluster_aligned_seq = alignment.get_read()
    alg_len = len(ref_aligned_seq)
    seq_to_align_len = len(seq_to_align)
    # Alignment reading parameters: from left to right for forward reads
    # from right to left for reverse reads
    if fwd:
        seq_pos, alg_pos = 0, 0
        inc_pos, stop_pos = 1, seq_to_align_len
    else:
        seq_pos, alg_pos = seq_to_align_len - 1, alg_len - 1
        inc_pos, stop_pos = -1, -1
    # Reading the alignment
    seq_aligned_list = [GAP for pos in range(alg_len)]
    ref_realigned_list = [GAP for pos in range(alg_len)]
    while seq_pos != stop_pos:
        if cluster_aligned_seq[alg_pos] != GAP:
            seq_aligned_list[alg_pos] = seq_to_align[seq_pos]
            seq_pos += inc_pos
        ref_realigned_list[alg_pos] = ref_aligned_seq[alg_pos]
        alg_pos += inc_pos
    # Trmming trailing gaps
    if fwd:
        seq_aligned = (''.join(seq_aligned_list)).rstrip(GAP)
        ref_realigned = (''.join(ref_realigned_list)).rstrip(GAP)
    else:
        seq_aligned = (''.join(seq_aligned_list)).lstrip(GAP)
        ref_realigned = (''.join(ref_realigned_list)).lstrip(GAP)
    return Alignment([ref_realigned, seq_aligned])


class Alignment:
    """
    Single alignment between a sequence  and an amplicon.

    :attribute: read (str): sequence of the read in the alignment (with gaps)
    :attribute: amplicon (str): amplicon sequence in the alignment (with gaps)
    :attribute: cigar (str): CIGAR string describing the alignment
    :attribute: normalized (bool): True if all variant blocks are normalized;
    initially set to None and updated when blocks are computed in variant_utils
    See variants_utils for the definition of normalized block.
    :attribute: score (float or None): alignment score
    """

    def __init__(self, alignment, cigar=None, normalized=None, alg_parameters=None):
        """
        Constructor

        :param: alignment: alignment object as returned by Bio.pairwise2
        :param: cigar (str or None): if not None, CIGAR string of the alignment
        otherwise it needs to be computed
        :param: normalized (bool or None): value for normalized attribute
        :param: alg_parameters (dict(str, float)): alignment parameters
        """
        self.amplicon = alignment[0]
        self.read = alignment[1]
        if cigar is not None:
            self.cigar = cigar
        else:
            self.__set_cigar()
        self.normalized = normalized
        if len(alignment) > 2:
            self.score = alignment[2]
        elif alg_parameters is not None:
            self.score = compute_alignment_score(self.amplicon, self.read, alg_parameters)
        else:
            self.score = None

    def get_score(self):
        """
        :return: score of the alignment (float or None)
        """
        return self.score

    def set_normalized(self, value):
        """
        Set the boolean value for self.normalized
        """
        self.normalized = value

    def check_normalized(self):
        """
        :return: bool or None: normalized value
        """
        return self.normalized

    def get_read(self):
        """
        :return: str: aligned read sequence
        """
        return self.read

    def get_amplicon(self):
        """
        :return: str: aligned amplicon sequence
        """
        return self.amplicon

    def get_cigar(self):
        """
        :return: the CIGAR string of the alignment
        """
        return self.cigar

    def __set_cigar(self):
        """
        Computes CIGAR string describing the alignment
        :return: str: CIGAR string
        """
        self.cigar = compute_alignment_cigar(self.amplicon, self.read)

    def amplicon_pos_to_read_pos(self, pos):
        """
        Computes the position in unaligned read aligned to the base
        in position pos in amplicon
        :param: pos (int): position in unaligned amplicon
        """
        return coord_to_aligned(self.amplicon, self.read, pos)

    def to_str(self, clusters, cluster_id, alg_id):
        """
        Export an alignment into a string in a FASTA-like format in a file
        >cluster_id:alignment_id cigar
        amplicon aligned sequence
        read aligned consensus sequence
        average quality of read sequence with gap positions shown as UNAMBIGUOUS
        alternate read sequence with gap and unambiguous positions shown as UNAMBIGUOUS
        alternate read quality with gap and unambiguous positions shown as UNAMBIGUOUS
        :param: clusters (ClusterReadsSet): clusters aligned
        :param: cluster_id, alg_id (str): ID of cluster and alignment
        :return: str
        """
        alg_pos = range(len(self.read))  # Alignment positions
        alg_pos_nogaps = [pos for pos in alg_pos if self.read[pos] != GAP]
        # shift: number of visited gaps in the aligned cluster sequence
        # shift_to_alg[pos] = corresp. alignment pos. for cluster pos. pos
        # shift_to_cluster[pos] = corresp.cluster pos. for alignment pos. pos
        shift_to_cluster, shift_to_alg, shift = {}, {}, 0
        for pos in alg_pos:
            shift += (1 if self.read[pos] == GAP else 0)
            shift_to_cluster[pos] = pos - shift
            shift_to_alg[pos - shift] = pos
        # Aligned consensus sequence quality
        cluster = clusters.get_cluster(cluster_id)
        cluster_qual = phred_to_str(cluster.get_qual())
        consensus_qual = [UNAMBIGUOUS for pos in alg_pos]
        for pos in alg_pos_nogaps:
            consensus_qual[pos] = cluster_qual[shift_to_cluster[pos]]
        # Alternate sequence (non-chosen ambiguous bases) and phred quality
        alternate_seq = [UNAMBIGUOUS for pos in alg_pos]
        alternate_qual = [UNAMBIGUOUS for pos in alg_pos]
        for (pos, base, qual) in cluster.get_ambiguous_bases():
            alternate_seq[shift_to_alg[pos]] = base
            alternate_qual[shift_to_alg[pos]] = phred_to_str([qual])
        # Output string, with trailing spaces removed to save disk space
        out_str = dedent(
            f"""
            >{cluster_id}:{alg_id} {self.cigar} {self.normalized} {self.score} {cluster.get_size()}
            {self.amplicon}
            {self.read}
            {''.join(consensus_qual)}
            {''.join(alternate_seq).rstrip()}
            {''.join(alternate_qual).rstrip()}"""
        )
        return out_str[1:]


class AlignmentsSet:
    """
    Set of alignments between a set of clusters and an amplicon, for clusters
    of size at least parameters[ALG_CLUSTER_MIN_SIZE].
    Alignments are kept only if they are amplicon-consistent: there is no indel
    or mismatch in the primers and the alignments covers the whole amplicon.

    :attribute: alignments: double dictionary of amplicon-consistent Alignment
    objects indexed by cluster ID then alignment ID for each cluster ID.
    """

    def __init__(self, log_parameters, alignments):
        """
        Constructor from a set of alignments for a set of clusters

        :param: log_parameters (dict(str, float or str)): log parameters
        if None, no generated log
        :param: alignments, double dictionary indexed by cluster ID then
        alignment ID of Alignment objects
        """
        self.alignments = alignments
        if log_parameters is not None:
            log = self.__init_log()
            for cluster_id, cluster_alignments in alignments.items():
                self.__update_log(log, alignments[cluster_id])
            self.__write_log(log, log_parameters)

    @staticmethod
    def __align(amplicon, cluster_seq, parameters):
        """
        Computes the alignments for a given cluster against its source amplicon

        :param: amplicon (Amplicon): amplicon object
        :param: cluster_seq (str): cluster consensus sequence
        :param: parameters (dict(float)): alignments parameters
        :return: dictionary indexed by alignment_id of Alignment objects
        """

        def check_amplicon_consistent(alg_amp_seq, alg_read_seq, amplicon):
            """
            Checks if aligned sequences form an amplicon-consistent alignment
            :param: amplicon (Amplicon): amplicon object
            :param: alg_amp_seq, alg_read_seq (str): aligned sequences
            :return: bool: True if alignment matches with primers and covers
            the whole amplicon sequence
            """
            test_read_primers = amplicon.check_primers(alg_read_seq)
            test_amp_primers = amplicon.check_primers(alg_amp_seq)
            test_full_length = alg_amp_seq.replace(GAP, '') == amplicon.get_seq()
            return test_read_primers and test_amp_primers and test_full_length

        # Alignment scoring scheme
        m, mm, go, ge = parameters[ALG_M], parameters[ALG_MM], parameters[ALG_GO], parameters[ALG_GE
                                                                                              ]
        # Alugning the sequences
        algs = pairwise2.align.globalms(amplicon.get_seq(), cluster_seq, m, mm, go, ge)
        # Extracting the alignments
        nb_max_alignments, nb_alignments, alignments = parameters[ALG_NB_MAX], 0, {}
        alg_id = 0
        while alg_id < len(algs) and nb_alignments < nb_max_alignments:
            # Excluding non-amplicon-consistent alignments
            if check_amplicon_consistent(algs[alg_id][0], algs[alg_id][1], amplicon):
                alignments[str(alg_id)] = Alignment(algs[alg_id], cigar=None)
                nb_alignments += 1
            alg_id += 1
        return alignments

    @classmethod
    def from_clusters(cls, log_parameters, amplicon, alg_parameters, clusters):
        """
        Creates an object from a set of clusters that are aligned against their
        source amplicon
        :param: amplicon (Amplicon): amplicon on which clusters are aligned
        :param: cluster (ClusterReadsSet): cluster: if not None, alignments are
        obtained from the clusters, otherwise they are imported from in_file
        :param: alg_parameters (dict(float)): alignments parameters
        """
        alignments = {}
        cluster_min_size = alg_parameters[ALG_CLUSTER_MIN_SIZE]
        amp_seq = amplicon.get_seq()
        # Loop on all clusters
        for cluster_id in clusters.get_clusters_id(to_sort=True):
            # Current cluster
            cluster = clusters.get_cluster(cluster_id)
            cluster_seq = cluster.get_consensus_seq()
            cluster_size = cluster.get_size()
            if cluster_size >= cluster_min_size and amp_seq != cluster_seq:
                # Different sequences: need to be aligned
                alignments[cluster_id] = cls.__align(amplicon, cluster_seq, alg_parameters)
            if cluster_size >= cluster_min_size and amp_seq == cluster_seq:
                # Skipping computing a trivial alignment between identical sequences
                alg, cigar = (amp_seq, cluster_seq), f"{CIGAR_MATCH}{len(amp_seq)}"
                alignment = Alignment(alg, cigar=cigar, alg_parameters=alg_parameters)
                alignments[cluster_id] = {'0': alignment}
        return cls(log_parameters, alignments)

    @classmethod
    def from_file(cls, log_parameters, in_file):
        """
        Creates an object from a file
        :param: log_parameters (dict(str, float or str)): log parameters
        if None, no generated log
        :param: in_file (str): path to alignments file
        """

        # Functions defined to process each line of imported file
        def __read_header(data_line):  # Alignment header
            data_split = data_line[1:].split()
            cluster_id, alg_id = data_split[0].split(':')
            cigar = data_split[1]
            str_norm = data_split[2]
            bool_norm = (None if str_norm == 'None' else bool(str_norm == 'True'))
            str_score = data_split[3]
            score = (None if str_score == 'None' else float(str_score))
            return cluster_id, alg_id, cigar, bool_norm, score

        # Lines 4, 5, 6 showing ambiguous bases are not used to read an alignment
        data_line_reader = {1: __read_header, 2: lambda x: x, 3: lambda x: x}

        alignments_data = open(in_file, 'r').readlines()
        line_id, data_line_values, alignments = 0, {}, defaultdict(dict)
        for data_line in alignments_data:
            line_id = (line_id % 6) + 1
            if data_line[0:4] != '#END' and line_id <= 3:
                data_line_values[line_id] = data_line_reader[line_id](data_line.rstrip())
            if line_id == 3:
                cluster_id, alg_id, cigar, bool_norm, score = data_line_values[1]
                amp_seq, read_seq = data_line_values[2], data_line_values[3]
                alg = [amp_seq, read_seq, score]
                alignment = Alignment(alignment=alg, cigar=cigar, normalized=bool_norm)
                alignments[cluster_id][alg_id] = alignment
        return cls(log_parameters, alignments)

    @staticmethod
    def __init_log():
        """
        0: nb alignments, 1: nb aligned clusters, 2: max number of alignments
        """
        return {0: 0, 1: 0, 2: 0}

    @staticmethod
    def __update_log(log, alignments):
        nb_algs = len(alignments.keys())
        log[0] += nb_algs
        if nb_algs > 0:
            log[1] += 1
        log[2] = max(nb_algs, log[2])

    @staticmethod
    def __write_log(log, parameters):
        sample_id, amplicon_id = parameters[LOG_SAMPLE], parameters[LOG_AMPLICON]
        log_pref = f"{sample_id}.{amplicon_id}"
        if log[1] > 0:
            avg_nb_algs = round(log[0] / log[1], 2)
            log_str = (
                f"{log_pref}\talignments:{log[0]} aligned_clusters:{log[1]} "
                f"max/average alignments:{log[2]} {avg_nb_algs}"
            )
            logger.info(log_str)
        else:
            logger.warning(f"{log_pref}\tno alignment")

    def get_alignments_id(self):
        """
        :return: dict(list(str)): list of all alignments ID indexed by
        clusters ID
        """
        return {c: list(self.alignments[c].keys()) for c in self.alignments.keys()}

    def get_alignment(self, cluster_id, alignment_id):
        """
        :return: Alignment: alignment of ID alignment_id for cluster cluster_id
        """
        return self.alignments[cluster_id][alignment_id]

    def export(self, out_file, clusters):
        """
        Print all amplicon-consistent alignments in a FASTA-like format in a file
        >cluster_id:alignment_id cigar
        amplicon aligned sequence
        read aligned consensus sequence
        average quality of read sequence with gap positions shown as UNAMBIGUOUS
        alternate read sequence with gap and unambiguous positions shown as UNAMBIGUOUS
        alternate read quality with gap and unambiguous positions shown as UNAMBIGUOUS

        :param: out_file (str): path to output file
        :param: clusters (ClusterReadsSet): clusters aligned
        """
        out, out_str = open(out_file, 'w'), ''
        for cluster_id in self.alignments.keys():
            for alg_id in self.alignments[cluster_id]:
                alignment = self.alignments[cluster_id][alg_id]
                out_str += f"\n{alignment.to_str(clusters, cluster_id, alg_id)}"
        out_str += '\n#END'
        out.write(out_str[1:])
        out.close()
