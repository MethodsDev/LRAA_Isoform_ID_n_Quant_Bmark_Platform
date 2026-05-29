
set -ex

if [ ! -f 'bmark.ok' ]; then

    $(git rev-parse --show-toplevel)/benchmarking/bmark_nb_runner.py \
        --analysisType ref_guided \
        --truth_gtf ../reference_data/GRCh38.gencode.v39.annotation.sorted.gtf \
        --truth_quant raw_prog_results/SGNex_A549_cDNAStranded_replicate5_run2.oarfish-v0.9.4.Oarfish.byAlignment.quant \
        --truth_reduced_gtf ../reference_data/GRCh38.gencode.v39.annotation_one_transcript_per_gene_with_exons_and_genes_ordered.gtf \
        --Venn_mode

    touch bmark.ok

fi



    
