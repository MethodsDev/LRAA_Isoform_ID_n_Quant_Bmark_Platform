#!/bin/bash

set -ex

if [ ! -f 'bmark.ok' ]; then

#    $(git rev-parse --show-toplevel)/benchmarking/bmark_nb_runner.py \
#        --analysisType ref_guided --truth_gtf ../reference_data/Araport11_GTF_genes_transposons.May2023.pigeon2_cleaned_with_gene.gtf --truth_quant ../reference_data/Arabidopsis_transcript_id_count.tsv --truth_reduced_gtf ../reference_data/Araport11_GTF_genes_transposons.May2023.pigeon2_cleaned_with_gene_1isoformPerGene_reduced.gtf

    $(git rev-parse --show-toplevel)/benchmarking/bmark_nb_runner.py \
        --analysisType ref_guided \
        --truth_gtf ../reference_data/Araport11_GTF_genes_transposons.May2023.pigeon2_cleaned_with_gene_expressed.gtf \
        --truth_quant ../reference_data/Arabidopsis_transcript_id_count.tsv \
        --truth_reduced_gtf ../reference_data/Araport11_GTF_genes_transposons.May2023.pigeon2_cleaned_with_gene_1isoformPerGene_expressed_kept.gtf
    
    touch bmark.ok

fi
