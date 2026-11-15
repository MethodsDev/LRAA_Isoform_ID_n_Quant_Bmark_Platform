#!/bin/bash
# Tool-neutral sets for the tool comparison. One library at a time; each step memory-capped.
RUN="systemd-run --user --scope -q -p MemoryMax=20G -p MemorySwapMax=0"
for s in BT474 K562 HG002 UHRR; do
  d=${s}_pool
  $RUN python3 make_pool_gtf.py $s && $RUN python3 enumerate_pairs.py $d ../out/pairs/$d.pairs.tsv.gz && $RUN python3 scan_reads.py $d 6 \
   && $RUN python3 score_sites.py $d ../out/scored/$d.sites.tsv.gz && $RUN python3 protein.py $d ../out/scored/$d.pairs.tsv.gz \
   && echo "$s done $(date +%H:%M)" || echo "$s FAILED $(date +%H:%M)"
done
echo POOL_DONE
