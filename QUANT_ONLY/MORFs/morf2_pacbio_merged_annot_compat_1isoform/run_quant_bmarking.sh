#!/bin/bash

set -ex


if [ ! -f 'bmark.ok' ]; then

    $(git rev-parse --show-toplevel)/benchmarking/bmark_nb_runner.py \
        --analysisType quant_only --truth_gtf reference_data/minigenome.UTRs_trimmed_1isoformref_expressed_pacbio.gtf --truth_quant reference_data/morf2_pacbio_merged_annot_compat_sorted_tn_counts.tsv

    touch bmark.ok

fi
