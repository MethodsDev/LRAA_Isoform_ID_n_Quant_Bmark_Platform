#!/bin/bash

set -ex

if [ ! -f 'bmark.ok' ]; then

    $(git rev-parse --show-toplevel)/benchmarking/bmark_nb_runner.py \
        --analysisType quant_only --truth_gtf ../reference_data/Araport11_GTF_genes_transposons.May2023.pigeon2_cleaned_with_gene.gtf --truth_quant ../reference_data/Arabidopsis_transcript_id_count.tsv

    touch bmark.ok

fi

