#!/bin/bash
# ORF prediction for the proteoform analysis (the paper's two proteoform supplementary figures) on one study's ref-guided,
# sc cluster-guided GTF, with pyTD (PyTransDecoder): cDNA from gffread, ORFs predicted on it and mapped to the
# genome through the GTF. Outputs go to orfs/<PBMC|FTD>/, where the Zenodo bundle puts the paper's ORFs.
#
#   run_orfs.sh <PBMC|FTD>
#
# Settings: orf_settings.conf. The GTF comes from the study submodule with its Zenodo bundle unpacked
# (fetch_data.sh). Needs gffread, pyTD and diamond on PATH (or PYTD_BIN).
set -euo pipefail
DS="${1:?usage: run_orfs.sh <PBMC|FTD>}"
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
case "$DS" in
  PBMC) GTF=$REPO/LRAA_paper-PBMCs_study/pbmcs_refGuided/pbmcs_refGuided_sc_cluster_guided/data/PBMCs_pbio_CG_refguided.withGeneSymbols.gtf
        P=PBMCs_pbio.LRAA.refGuided.scg ;;
  FTD)  GTF=$REPO/LRAA_paper-FTD_study/als_refGuided/als_refGuided_sc_cluster_guided/data/LP00714.refguided.withGeneSymbols.gtf
        P=LP00714.refGuided.scg ;;
  *) echo "unknown dataset $DS" >&2; exit 2 ;;
esac
[ -s "$GTF" ] || { echo "missing $GTF: unpack the study's Zenodo bundle first" >&2; exit 3; }
# shellcheck source=/dev/null
source "$HERE/orf_settings.conf"
OUT=$HERE/orfs/$DS; mkdir -p "$OUT"; cd "$OUT"
PYTD_VERSION=$("$PYTD_BIN" --version 2>&1 | awk '{print $NF}')
SETTINGS="ORF_TOOL=pyTD PYTD_VERSION=$PYTD_VERSION MIN_PROT_LEN=$MIN_PROT_LEN SINGLE_BEST_ONLY=$SINGLE_BEST_ONLY MAX_QVALUE=$MAX_QVALUE HOMOLOGY_TOOL=$HOMOLOGY_TOOL HOMOLOGY_DB=$(basename "$HOMOLOGY_DB") GENOME=$(basename "$GENOME")"
ARGS=(-t "$P.cdna.fasta" --structures "$GTF" -O . -m "$MIN_PROT_LEN" --max-qvalue "$MAX_QVALUE"
      --blast-search-pep "$HOMOLOGY_DB" --blast-tool "$HOMOLOGY_TOOL" --blast-threads "$THREADS"
      --cpu "$THREADS")
[ "$SINGLE_BEST_ONLY" = 1 ] || ARGS+=(--retain-multiple-orfs)
echo "$DS: ORF prediction on $(basename "$GTF"), $SETTINGS, $(date '+%F %T')"
# transcripts named by transcript_id; LRAA writes every strand as + or -
gffread -w "$P.cdna.fasta" -g "$GENOME" "$GTF"
"$PYTD_BIN" "${ARGS[@]}" > run_orfs.log 2>&1
for f in "$P.cdna.fasta.transdecoder.pep" "$P.cdna.fasta.transdecoder.genome.gff3"; do
    [ -s "$f" ] || { echo "missing output $f, see $OUT/run_orfs.log" >&2; exit 1; }
done
printf "LRAA_VERSION=v0.51.0\nGTF=%s\nSETTINGS=%s\nDATE=%s\n" "$(basename "$GTF")" "$SETTINGS" "$(date '+%F %T')" > ORFS_VERSION.txt
echo "$DS ORFs done"
