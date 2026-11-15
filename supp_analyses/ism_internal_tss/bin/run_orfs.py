#!/usr/bin/env python3
"""ORF route for L models with no GENCODE chain match, for S sites with read pile-up and cap evidence.
Writes L models to a GTF, extracts their cDNA with gffread and predicts ORFs with pyTD (PyTransDecoder), using the
proteoform analysis's shared settings (LRAA_proteoform_analyses/orf_settings.conf: pyTD v0.4.0, >= 50 aa, single best
ORF, diamond vs GENCODE v47 proteins at e-value 1e-5, target-decoy q <= 0.05 for ORFs without a homology hit, start
refinement on), plus complete ORFs only (this study needs L's start codon). ORFs are placed on the genome through the
GTF (--structures). Writes out/orf/orfs.tsv: chrom strand l_chain cds. usage: run_orfs.py CANDIDATES_TSV"""
import sys, subprocess
import pandas as pd
from common import *
here = os.path.dirname(os.path.abspath(__file__)); W = os.path.abspath(os.path.join(here, "..", "out", "orf")); os.makedirs(W, exist_ok=True)
d = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)
d = d[(d.cds_route == "L_unannotated") & d.ev_pile & d.ev_cap & d["class"].isin(["strict_ISM", "TSS_in_L_intron"])]
L = d.sort_values("l_up_junc_reads", ascending=False).drop_duplicates(["chrom", "strand", "l_chain"])
print(len(d), "pairs", len(L), "L chains", flush=True)
with open(f"{W}/L.gtf", "w") as o:
    for i, r in enumerate(L.itertuples()):
        c = sorted(parse_chain(r.l_chain))
        lo, hi = (r.l_tss, r.s_end3) if r.strand == "+" else (r.s_end3, r.l_tss)
        ex = [(lo, c[0][0] - 1)] + [(c[j][1] + 1, c[j + 1][0] - 1) for j in range(len(c) - 1)] + [(c[-1][1] + 1, hi)]
        o.write(f'{r.chrom}\tLRAA\ttranscript\t{lo}\t{hi}\t.\t{r.strand}\t.\tgene_id "g{i}"; transcript_id "L{i}";\n')
        for a, b in ex:
            o.write(f'{r.chrom}\tLRAA\texon\t{a}\t{b}\t.\t{r.strand}\t.\tgene_id "g{i}"; transcript_id "L{i}";\n')
L[["chrom", "strand", "l_chain"]].assign(tx=[f"L{i}" for i in range(len(L))]).to_csv(f"{W}/L.index.tsv", sep="\t", index=False)
# shared ORF settings (bash KEY=VALUE file; $CTAT_GENOME_LIB and ~ expanded as run_orfs.sh does)
os.environ.setdefault("CTAT_GENOME_LIB", f"{CTAT}/ctat_genome_lib_build_dir")
conf = os.path.join(REPO, "LRAA_proteoform_analyses", "orf_settings.conf")
cfg = {}
for line in open(conf):
    m = re.match(r"^([A-Z_]+)=(\S+)", line)
    if m: cfg[m.group(1)] = os.path.expandvars(os.path.expanduser(m.group(2)))
for f in os.listdir(W):                      # outputs of an earlier run, so nothing stale is read
    if f.startswith("L.") and f not in ("L.gtf", "L.index.tsv") or f == "__chckpts":
        subprocess.check_call(["rm", "-rf", os.path.join(W, f)])
subprocess.check_call(["gffread", "-w", "L.cdna.fasta", "-g", SC_GENOME, "L.gtf"], cwd=W)
args = [cfg["PYTD_BIN"], "-t", "L.cdna.fasta", "--structures", "L.gtf", "-O", ".", "-m", cfg["MIN_PROT_LEN"], "--max-qvalue", cfg["MAX_QVALUE"],
        "--blast-search-pep", cfg["HOMOLOGY_DB"], "--blast-tool", cfg["HOMOLOGY_TOOL"], "--blast-threads", cfg["THREADS"], "--cpu", cfg["THREADS"],
        "--complete-orfs-only"]
if cfg["SINGLE_BEST_ONLY"] != "1": args.append("--retain-multiple-orfs")
ver = subprocess.run([cfg["PYTD_BIN"], "--version"], capture_output=True, text=True).stdout.split()[-1]
print("pyTD", ver, " ".join(args), flush=True)
subprocess.check_call(args, cwd=W)
open(f"{W}/ORFS_VERSION.txt", "w").write(f"ORF_TOOL=pyTD\nPYTD_VERSION={ver}\nSETTINGS={conf}\nCMD={' '.join(args)}\n")
idx = pd.read_csv(f"{W}/L.index.tsv", sep="\t").set_index("tx")
cds = {}
for line in open(f"{W}/L.cdna.fasta.transdecoder.genome.gff3"):
    p = line.rstrip().split("\t")
    if len(p) > 8 and p[2] == "CDS":
        tx = re.search(r"Parent=(?:cds\.)?(L\d+)", p[8]).group(1)
        cds.setdefault(tx, []).append(f"{p[3]}-{p[4]}")
with open(f"{W}/orfs.tsv", "w") as o:
    o.write("chrom\tstrand\tl_chain\tcds\n")
    for t, v in cds.items():
        r = idx.loc[t]
        o.write(f"{r.chrom}\t{r.strand}\t{r.l_chain}\t{','.join(sorted(v, key=lambda x: int(x.split('-')[0])))}\n")
print("ORFs", len(cds), "of", len(L))
