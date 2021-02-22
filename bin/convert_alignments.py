#!/usr/bin/env python3
"""
Converting alignments text files into BAM files
"""

# Standard imports
import os
from textwrap import dedent
import sys

# Third-party imports
from Bio import Seq, SeqIO
from Bio.Sequencing.Applications import SamtoolsViewCommandline

# Local imports
from alignments_utils import (
    align_single_read,
    AlignmentsSet,
    GAP,
    compute_alignment_cigar,
    compute_clipped_alignment_cigar,
    compute_alignment_MD,
    compute_alignment_NM,
)
from bin.amplicons_utils import AmpliconsManifest
from bin.clusters_utils import ClusterReadsSet
from bin.fastq_utils import phred_to_str

# Auxiliary functions related to SAM files

SAM_HEADER_SQ = dedent(
    f"""
    @SQ\tSN:chrM\tLN:16571
    @SQ\tSN:chr1\tLN:249250621
    @SQ\tSN:chr2\tLN:243199373
    @SQ\tSN:chr3\tLN:198022430
    @SQ\tSN:chr4\tLN:191154276
    @SQ\tSN:chr5\tLN:180915260
    @SQ\tSN:chr6\tLN:171115067
    @SQ\tSN:chr7\tLN:159138663
    @SQ\tSN:chr8\tLN:146364022
    @SQ\tSN:chr9\tLN:141213431
    @SQ\tSN:chr10\tLN:135534747
    @SQ\tSN:chr11\tLN:135006516
    @SQ\tSN:chr12\tLN:133851895
    @SQ\tSN:chr13\tLN:115169878
    @SQ\tSN:chr14\tLN:107349540
    @SQ\tSN:chr15\tLN:102531392
    @SQ\tSN:chr16\tLN:90354753
    @SQ\tSN:chr17\tLN:81195210
    @SQ\tSN:chr18\tLN:78077248
    @SQ\tSN:chr19\tLN:59128983
    @SQ\tSN:chr20\tLN:63025520
    @SQ\tSN:chr21\tLN:48129895
    @SQ\tSN:chr22\tLN:51304566
    @SQ\tSN:chrX\tLN:155270560
    @SQ\tSN:chrY\tLN:59373566"""
)[1:]

SAM_VERSION = '1.6'

def SAM_header():
    sam_version = f"VN:{SAM_VERSION}"
    sam_sorting = f"SO:unsorted"
    sam_grouping = f"GO:none"
    HD_info = '\t'.join([sam_version, sam_sorting, sam_grouping])
    HD = f"@HD\t{HD_info}"
    SQ = SAM_HEADER_SQ
    PG = f"@PG\tID:Canexia indels caller"
    return '\n'.join([HD, SQ, PG])


def compute_flag(fwd=True, primary=True, unmapped=False):
    """
    Compute the SAM FLAG for an alignment
    :param: fwd (bool): True if forward read
    :param: primary (bool): True if first encountered alignment for this read
    :param: unmapped (bool): True if read unmapped
    :return: int

    unmapped reads: 77/141
      read paired (0x1)
      read unmapped (0x4)
      mate unmapped (0x8)
      first in pair (0x40) --> forward read
      second in pair (0x80) --> reverse read

    fwd mapped reads: 99/355
      read paired (0x1)
      read mapped in proper pair (0x2)
      mate reverse strand (0x20)
      first in pair (0x40)
      not primary alignment (0x100) --> for all but the first alignment

    rev mapped: 147/403
      read paired (0x1)
      read mapped in proper pair (0x2)
      read reverse strand (0x10)
      second in pair (0x80)
      not primary alignment (0x100) --> for all but the first alignment
    """
    if unmapped:
        return (77 if fwd else 141)
    if fwd:
        return (99 if primary else 355)
    return (147 if primary else 403)



def alignment_to_sam(read_id, alignment, amplicon, fwd=True, primary=True):
    """
    Computes the SAM string for an aligned read
    :param: read_id (str): read ID
    :param: alignment (str, str, str, str, float): read aligned sequence,
    read Phred quality string, reference aligned sequence, CIGAR string, score
    of the alignment
    :param: amplicon (Amplicon): amplicon object where the read is aligned
    :param: fwd (bool): True if forward read, false if reverse read
    :param: primary (bool): True if first encountered alignment for the read
    :return: str: SAM string
    """
    # Reading from alignment the elements required to realign a read
    read_aligned_seq = alignment[0]
    ref_realigned_seq = alignment[1]
    read_qual, read_cigar = alignment[1], alignment[3]
    score, NM, MD = alignment[4], alignment[5], alignment[6]
    # Populating SAM fields
    QNAME = read_id
    FLAG = compute_flag(fwd=fwd, primary=primary, unmapped=False)
    RNAME = amplicon.get_chr()
    POS = amplicon.get_start() + 1 # +1 because pos is 1-based in SAM
    MAPQ = 255
    CIGAR = read_cigar
    RNEXT = '='
    PNEXT = POS
    amplicon_len = amplicon.get_len()
    TLEN = (amplicon_len if fwd else -amplicon_len)
    SEQ = read_aligned_seq.replace(GAP, '')
    QUAL = read_qual
    AS = f"AS:i:{score}"
    NM = f"NM:i:{NM}"
    MD = f"MD:Z:{MD}"
    SAM_fields = [
        QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL,
        NM, MD, AS
    ]
    SAM_str = '\t'.join([str(x) for x in SAM_fields])
    return  SAM_str


# Handling mapped merged reads
# Main function: process_read generates the SAM entry for a pair of forward and
# reverse read originating from a merged read

def get_reads_info(read_id, fwd_fq_index, rev_fq_index):
    """
    Given a read ID and the index of the foward reads and reverse reads
    returns the forward read sequence and the reverse complement of the reverse
    read
    :param: read_id (str): FASTQ read ID
    :param: fwd_fq_index, rev_fq_index: indexes for the forward and
    reverse reads files
    :return: str, str, str, str: forward read sequence and quality, reverse
    complement of the reverse read, reverse of the quality of the reverse read
    """
    fwd_record = fwd_fq_index[read_id]
    fwd_seq = fwd_record.seq
    fwd_qual = phred_to_str(fwd_record.letter_annotations["phred_quality"])
    rev_record = rev_fq_index[read_id]
    rev_seq = Seq.reverse_complement(rev_record.seq)
    rev_qual_int = reversed(rev_record.letter_annotations["phred_quality"])
    rev_qual = phred_to_str(rev_qual_int)
    return fwd_seq, fwd_qual, rev_seq, rev_qual

def trim_read_seq(cluster_seq, read_seq, fwd=True):
    """
    Separate a read into the part belonging to the merged read used to define
    the cluster consensus sequence cluster_seq is trimmed.
    If fwd is True, assumes read_seq is a forward read and the trimmed part is
    a suffix.
    If fwd is False, assumes read_seq is the reverse complement of a reverse
    read and the trimmed part is a prefix.
    :param: cluster_seq (str): cluster consensus sequence
    :param: read_seq (str): read sequence
    :param: fwd (bool): if True, read_sew is assumed to be a forward read
    otherwise a reverse read.
    :return: str, str: prefix, suffix of the trimmed read.
    """
    cluster_len, read_len = len(cluster_seq), len(read_seq)
    common_len = min(cluster_len, read_len)
    if fwd:
        read_prefix_seq = read_seq[:common_len]
        read_suffix_seq = read_seq[common_len:]
    else:
        read_prefix_seq = read_seq[:- common_len]
        read_suffix_seq = read_seq[- common_len:]
    return read_prefix_seq, read_suffix_seq

def process_single_read(alignment, cluster_seq, read_seq, read_qual, fwd=True):
    """ Given a read, an alignment and a cluster consensus sequence, realign
    the read sequence according to the alignment and returns the aligned
    sequences and the alignment cigar string.
    :param: alignment (Alignment): cluster alignment
    :param: cluster_seq (str): cluster consensus sequence
    :param: read_seq (str): full read sequence
    :param: read_qual (str): Phred quality string
    :param: fwd (bool): True if forward read
    :return: str, str, str, str: read aligned sequence, read quality over the
    aligned part, reference realigned sequence, cigar string
    """
    read_prefix_seq, read_suffix_seq = trim_read_seq(
        cluster_seq, read_seq, fwd=fwd
    )
    if fwd:
        read_to_align_seq, read_clipped_seq = read_prefix_seq, read_suffix_seq
        trimmed_read_qual = read_qual[:len(read_to_align_seq)]
    else:
        read_to_align_seq, read_clipped_seq = read_suffix_seq, read_prefix_seq
        trimmed_read_qual = read_qual[-len(read_to_align_seq):]
    read_alignment = align_single_read(alignment, read_to_align_seq, fwd=fwd)
    read_aligned_seq = read_alignment.get_read()
    ref_realigned_seq = read_alignment.get_amplicon()
    read_cigar = compute_clipped_alignment_cigar(
        ref_realigned_seq, read_aligned_seq, read_clipped_seq, fwd=fwd
    )
    score = read_alignment.get_score()
    NM = compute_alignment_NM(ref_realigned_seq, read_aligned_seq)
    MD = compute_alignment_MD(ref_realigned_seq, read_aligned_seq)
    return (
        read_aligned_seq, trimmed_read_qual, ref_realigned_seq,
        read_cigar, score, NM, MD
    )

def process_merged_read(
    alignment, cluster_seq, read_id, fwd_fq_index, rev_fq_index, primary=False
):
    """ Given a read_id, an alignment and a cluster consensus sequence, realign
    the forward and reverse reads sequences according to the alignment and
    returns the corresponding SAM strings.
    :param: alignment (Alignment): cluster alignment
    :param: cluster_seq (str): cluster consensus sequence
    :param: read_id (str): ID of the merged read
    :param:  fwd_fq_index, rev_fq_index: indexes for the forward and
    reverse reads files
    :return: str: SAM string, one per line, for the alignments of the forward
    and reverse reads constrained by the alignment of the merged read
    """
    fwd_seq, fwd_qual, rev_seq, rev_qual = get_reads_info(
        read_id, fwd_fq_index, rev_fq_index
    )
    fwd_alignment = process_single_read(
        alignment, cluster_seq, fwd_seq, fwd_qual, fwd=True
    )
    rev_alignment = process_single_read(
        alignment, cluster_seq, rev_seq, rev_qual, fwd=False
    )
    SAM_str = '\n' + alignment_to_sam(
        read_id, fwd_alignment, amplicon,
        fwd=True, primary=primary
    )
    SAM_str += '\n' + alignment_to_sam(
        read_id, rev_alignment, amplicon,
        fwd=False, primary=primary
    )
    return SAM_str

# Handling unmapped (unmerged) reads
# Main function: process_unmapped_read generates the SAM entry for the a pair of
# forward and reverse reads that were not mapped from a merged read.

def get_unmerged_reads_id(unmapped_reads_file):
    """
    Reads the unmerged reads file and returns a list of read ID
    :param: unmapped_reads_file (str): path to file recording all unmapped reads
    ID (one ID per pair per line, first field)
    :return: list(str): list of read ID
    :TODO: integrate this to reads_utils as a non-class function
    """
    reads_id_data = open(unmapped_reads_file, 'r').readlines()
    return [x.rstrip().split('\t')[0] for x in reads_id_data]

def unmapped_to_SAM(read_id, read_seq, read_qual, fwd=True):
    """
    Computes the SAM string for an aligned read
    :param: read_id (str): read ID
    :param: read_seq (str): read sequence
    :param: read_qual (str): read Phred quality string
    :param: fwd (bool): True if forward read, false if reverse read
    :return: str: SAM string
    """
    QNAME = read_id
    FLAG = compute_flag(fwd=fwd, primary=False, unmapped=True)
    RNAME = amplicon.get_chr()
    POS = 0
    MAPQ = 255
    CIGAR = '*'
    RNEXT = '*'
    PNEXT = 0
    TLEN = 0
    SEQ = read_seq
    QUAL = read_qual
    SAM_fields = [
        QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL
    ]
    return  '\t'.join([str(x) for x in SAM_fields])

def process_unmapped_reads(unmapped_reads_file, fwd_fq_index, rev_fq_index):
    """
    Compute the SAM string for all unmapped reads
    :param: unmapped_reads_file (str): path to file recording all unmapped reads
    ID (one ID per pair)
    :param:  fwd_fq_index, rev_fq_index: indexes for the forward and
    reverse reads files
    :return: str: SAM string, one per line, for each read
    """
    reads_id = get_unmerged_reads_id(unmapped_reads_file)
    SAM_str = ''
    for read_id in reads_id:
        fwd_seq, fwd_qual, rev_seq, rev_qual = get_reads_info(
            read_id, fwd_fq_index, rev_fq_index
        )
        SAM_str += '\n'
        SAM_str += f"{unmapped_to_SAM(read_id, fwd_seq, fwd_qual, fwd=True)}"
        SAM_str += '\n'
        SAM_str += f"{unmapped_to_SAM(read_id, rev_seq, rev_qual, fwd=False)}"
    return SAM_str

# Main function, currently aimed at testing the code.

if __name__ == "__main__":
    # Parameters
    # Directory containing alignments and clusters text files, and fastq files
    DATA_DIR = sys.argv[1]
    # List of samples to test the code on
    SAMPLES_FILE = os.path.join(DATA_DIR, 'samples.txt')
    # List of amplicons to test the code on
    AMPLICONS_FILE = os.path.join(DATA_DIR, 'amplicons.txt')
    # Amplicons manifest
    MANIFEST_FILE_NAME = 'CG001v4.0_Amplicon_Manifest_Panel4.0.3_20181101.tsv'
    MANIFEST_FILE_PATH = os.path.join(DATA_DIR, MANIFEST_FILE_NAME)

    # Reading input data
    AMPLICON_MANIFEST = AmpliconsManifest(MANIFEST_FILE_PATH )
    SAMPLES_ID = [x.rstrip()  for x in open(SAMPLES_FILE, 'r').readlines()]
    AMPLICONS_ID = [x.rstrip()  for x in open(AMPLICONS_FILE, 'r').readlines()]

    # Processing samples
    for sample_id in SAMPLES_ID:
        # Reading annotated FASTQ files
        FQ_FWD_FILE_NAME = f"{sample_id}_L001_R1_001_annotated.fq"
        FQ_FWD_FILE = os.path.join(DATA_DIR,FQ_FWD_FILE_NAME)
        FQ_FWD = SeqIO.index(FQ_FWD_FILE, "fastq")
        FQ_REV_FILE_NAME = f"{sample_id}_L001_R2_001_annotated.fq"
        FQ_REV_FILE = os.path.join(DATA_DIR, FQ_REV_FILE_NAME)
        FQ_REV = SeqIO.index(FQ_REV_FILE, "fastq")
        # Output SAM file
        SAM_out_file_name = f"{sample_id}_L001_001_annotated.sam"
        SAM_out_file_path = os.path.join(DATA_DIR, SAM_out_file_name)
        SAM_out = open(SAM_out_file_path, 'w')
        SAM_out.write(SAM_header())
        # Reading clusters and alignments
        for amplicon_id in AMPLICONS_ID:
            SAM_str = ''
            amplicon = AMPLICON_MANIFEST.get_amplicon(amplicon_id)
            # Reading clusters
            clusters_file_name = f"{sample_id}_{amplicon_id}_clusters.txt"
            clusters_file_path = os.path.join(DATA_DIR, clusters_file_name)
            clusters = ClusterReadsSet.from_file(None, clusters_file_path)
            clusters_id = clusters.get_clusters_id()
            # Reading alignments
            alignments_file_name = f"{sample_id}_{amplicon_id}_alignments.txt"
            alignments_file_path = os.path.join(DATA_DIR, alignments_file_name)
            alignments = AlignmentsSet.from_file(None, alignments_file_path)
            alignments_id = alignments.get_alignments_id()
            # Loop on alignments in the alignments text file
            for cluster_id in clusters_id:
                cluster = clusters.get_cluster(cluster_id)
                cluster_seq = cluster.get_consensus_seq()
                # Reads id of all reads in the cluster
                reads_id = cluster.get_reads_id()
                # Number of read alignments
                reads_alg_count = {read_id: 0 for read_id in reads_id}
                for alg_id in alignments_id[cluster_id]:
                    alignment = alignments.get_alignment(cluster_id, alg_id)
                    for read_id in reads_id:
                        primary = (reads_alg_count[read_id] == 0)
                        SAM_str += process_merged_read(
                            alignment, cluster_seq, read_id,
                            FQ_FWD, FQ_REV, primary=primary
                        )
                        reads_alg_count[read_id] += 1
            # Processing unmapped reads
            unmapped_reads_file_name = (
               f"{sample_id}_{amplicon_id}_unmerged_reads.txt"
            )
            unmapped_reads_file_path = os.path.join(
               DATA_DIR, unmapped_reads_file_name
            )
            SAM_str += process_unmapped_reads(
               unmapped_reads_file_path, FQ_FWD, FQ_REV
            )
            SAM_out.write(SAM_str)
        SAM_out.close()
        # Conversion of the temporary SAM file into a BAM file
        BAM_out_file_path = f"{SAM_out_file_path}.bam"
        samtools_view_cmd = SamtoolsViewCommandline(
            input_file=SAM_out_file_path,
            o=BAM_out_file_path
        )
        samtools_view_cmd.__call__()
