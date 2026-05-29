#!/bin/bash

set -ex

if [ ! -f 'bmark.ok' ]; then

    $(git rev-parse --show-toplevel)/benchmarking/bmark_nb_runner.py \
        --analysisType quant_only --truth_gtf ../reference_data/gencode.vM32.annotation.gtf --truth_quant ../reference_data/Mouse_isoseqsim_ground_truth_transcript_id_count.tsv

    touch bmark.ok

fi

