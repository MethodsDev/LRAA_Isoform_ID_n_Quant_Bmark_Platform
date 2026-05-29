set -ex

if [ ! -f 'bmark.ok' ]; then

    $(git rev-parse --show-toplevel)/benchmarking/bmark_nb_runner.py \
        --truth_gtf ../reference_data/GRCh38.gencode.v39.annotation.sorted.gtf \
        --analysisType quant_only_no_truthset

    touch bmark.ok

fi
