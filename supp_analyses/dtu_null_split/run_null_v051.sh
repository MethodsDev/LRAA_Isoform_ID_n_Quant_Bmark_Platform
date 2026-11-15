#!/bin/bash
# Null split on LRAA v0.51.0 (PBMC, ref-guided sc cluster-guided models, pooled final quant): the paper's
# feature-usage DEXSeq pipeline (settings and inputs of the PBMC study's __Splice_Pattern_Usage_Analysis), run over
# random half-clusters of clusters 0, 1, 4 and 10 (clusters.null_s<seed>.tsv: half b = 100 + cluster; other cells
# left out). Every call between two halves of one cluster is false.
# Inputs: the PBMC study submodule with its Zenodo bundle unpacked (fetch_data.sh), an LRAA checkout at
# LRAA_v0.51.0 (LRAA_HOME), and the GRCh38 genome FASTA the LRAA runs used (GENOME_FA).
# usage: run_null_v051.sh <seed> [cores]
set -euo pipefail
SEED=$1; CPU=${2:-8}
HERE=$(cd "$(dirname "$0")" && pwd)
E=$(cd "$HERE/../../LRAA_paper-PBMCs_study" && pwd)
UTIL=${LRAA_HOME:-$HOME/LongReadAlignmentAssembler}/util/sc/site_usage
GENOME=${GENOME_FA:?set GENOME_FA to the GRCh38 genome FASTA}
CG=$E/pbmcs_refGuided/pbmcs_refGuided_sc_cluster_guided
DATA=$CG/data
TX_MAP=$DATA/PBMCs_pbio.gene_transcript_splicehashcode.withGeneSymbols.tsv
QUANT_TAR=$DATA/PBMCs_pbio_CG_refguided.LRAA.final.cluster_quant.EXPRs.tar.gz
SP_MATRIX="$CG/scg_splice_collapsed/PBMCs_pbio.LRAA.sc^splice_pattern-sparseM"
SP_GTF=$CG/scg_splice_collapsed/LRAA.pbmcs.sc-cluster-guided.refGuided.sp_collapsed.gtf
ISO_MATRIX="$CG/diff_iso_usage_for_alt_termini_analysis/PBMCs_pbio.LRAA.sc^isoform-sparseM"
SITES=$E/__Site_Usage_Analysis/PBMCs.sites.tsv
POLYASITE=$E/__Site_Usage_Analysis/ext/atlas.clusters.2.0.GRCh38.96.bed.gz
CF=$HERE/clusters.null_s$SEED.tsv
OUT=$HERE/work/s$SEED; mkdir -p $OUT; cd $OUT
P=PBMCs; DEX=$P.dexseq
python3 $UTIL/feature_counts_from_sparse_matrix.py --mode splice_pattern --sparseM_dir "$SP_MATRIX" \
    --gene_transcript_splicehash $TX_MAP --cell_clusters $CF --sp_gtf $SP_GTF --output_prefix $P.SplicePattern > counts.SplicePattern.log 2>&1
python3 $UTIL/feature_counts_from_sparse_matrix.py --mode isoform_termini --sparseM_dir "$ISO_MATRIX" \
    --gene_transcript_splicehash $TX_MAP --cell_clusters $CF --sites $SITES --output_prefix $P.IsoformTermini > counts.IsoformTermini.log 2>&1
for K in SplicePattern IsoformTermini; do
  Rscript $UTIL/site_usage_dexseq.R --counts_prefix $P.$K --sites $P.$K.features.tsv --kind $K --cores $CPU \
      --output_prefix $DEX 2>&1 | grep -v "converting counts to integer mode" > dexseq.$K.log
  test -s $DEX.$K.pairwise.tsv.gz || { echo "no pairwise for $K"; continue; }
  python3 $UTIL/annotate_feature_usage_events.py --dexseq_prefix $DEX --kind $K --features $P.$K.features.tsv \
      --library $P.$K.cluster_library.tsv --cluster_quant_tar $QUANT_TAR --sites $SITES \
      --polyasite_atlas $POLYASITE --genome_fa $GENOME --output $DEX.$K.events.tsv > events.$K.log 2>&1
done
echo "$OUT DONE"
