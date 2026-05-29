
set -ex

if [ ! -f 'bmark.ok' ]; then 

    $(git rev-parse --show-toplevel)/benchmarking/bmark_nb_runner.py \
    --analysisType ref_free \
    --truth_gtf ../reference_data/GRCh38.gencode.v39.annotation.sorted.gtf \
    --truth_quant raw_prog_results/CL_BT474_E0_human.oarfish-v0.9.4.Oarfish.byAlignment.quant  \
    --Venn_mode

    touch bmark.ok

fi
