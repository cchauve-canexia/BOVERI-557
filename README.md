# BOVERI-557
Export alignments in BAM format

The indels pipeline generates alignments in text format.
In order to be visualized with IGV they should be also exported in BAM format.

Specification/Questions:
- for co-optimal alignments, the primary one is chosen arbitrarily (first alignment),
- one single BAM file per sample,
- additional file linking indel calls to alignments ID,
- alignments of original paired-end reads, not merged reads.
