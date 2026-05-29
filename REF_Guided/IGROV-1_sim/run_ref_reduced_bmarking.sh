#!/bin/bash

set -ex


$(git rev-parse --show-toplevel)/benchmarking/bmark_nb_runner.py \
    --analysisType ref_guided --truth_gtf reference_data/benchmark_full_annotations.wGene.gtf --truth_quant reference_data/benchmark_transcript_expression.tab --truth_reduced_gtf reference_data/benchmark_downsampled_30.wGene.gtf
