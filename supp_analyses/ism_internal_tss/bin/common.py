"""Dataset table and shared helpers for the ISM internal-TSS study. All paths relative to the repo root."""
import os, re, gzip, hashlib, collections
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
LRAA_HOME = os.path.join(REPO, "__LRAA_local_runs/v051_work/LRAA_v051")  # LRAA v0.51.0 code (de5e127b, detached worktree): read rules, ORF wrapper
CTAT = os.path.expanduser("~/CTAT_GENOME_LIBS/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play")
SC_GENOME = f"{CTAT}/ctat_genome_lib_build_dir/ref_genome.fa"
GENCODE47 = f"{CTAT}/GTEx_genome_resources/human.refdata-gex-GRCh38-GENCODE47.gtf"
BULK_GENOME = f"{REPO}/DENOVO_ID/CellLines/reference_data/GRCh38_no_alt.fa"
RDNA = f"{REPO}/__shared_refs/rdna_mask.GRCh38.bed"
EA = f"{REPO}/__Benchmark_extras/isoform_end_accuracy/data"
CAGE, REFTSS = f"{EA}/cage.interval.bed", f"{EA}/reftss.interval.bed"
P = f"{REPO}/PBMCs_kinnex/LRAA_PBMCs_eval"
A = f"{REPO}/ALS_LP00714/ALS_LRAA_eval"
BI = f"{REPO}/__LRAA_local_runs/bmark_v0.35.0_big/inputs"
TOOL = "LRAA-v0.51.0-de5e127"
# S and L share a 3' end when their last exons end within this many nt (user, 2026-11-11). Without it most strict ISM
# pairs in the single-cell sets are S models that run past L's 3' end, which LRAA never treats as contained.
END3_TOL = 50

DATASETS = {
 "PBMC": dict(kind="sc", bam=f"{REPO}/PBMCs_kinnex/PBMCs_pbio.aligned.sorted.bam", genome=SC_GENOME, baseline=0,
   clusters=f"{P}/pbmcs_refGuided/data/RefQuantOnly.cell_clusters.tsv", site_prefix=f"{P}/__Site_Usage_Analysis/PBMCs",
   gtfs={("denovo", "basic"): f"{P}/pbmcs_denovo/data/PBMCs_pbio_CG_denovo.LRAA.ref-free.gtf",
         ("denovo", "scg"): f"{P}/pbmcs_denovo/pbmcs_denovo_sc_cluster_guided/data/PBMCs_pbio_CG_denovo.withGeneSymbols.gtf",
         ("refGuided", "basic"): f"{P}/pbmcs_refGuided/data/PBMCs_pbio_CG_refguided.LRAA.ref-guided.gtf",
         ("refGuided", "scg"): f"{P}/pbmcs_refGuided/pbmcs_refGuided_sc_cluster_guided/data/PBMCs_pbio_CG_refguided.withGeneSymbols.gtf"}),
 "ALS": dict(kind="sc", bam=f"{REPO}/ALS_LP00714/LP00714.aligned.sorted.bam", genome=SC_GENOME, baseline=0,
   clusters=f"{A}/als_refGuided/data/RefQuantOnly.LP00714.ALS.cell_clusters.tsv", site_prefix=f"{A}/__Site_Usage_Analysis/LP00714",
   gtfs={("denovo", "basic"): f"{A}/als_denovo/data/LP00714.denovo.LRAA.ref-free.gtf",
         ("denovo", "scg"): f"{A}/als_denovo/als_denovo_sc_cluster_guided/data/LP00714.denovo.withGeneSymbols.gtf",
         ("refGuided", "basic"): f"{A}/als_refGuided/data/LP00714.refguided.LRAA.ref-guided.gtf",
         ("refGuided", "scg"): f"{A}/als_refGuided/als_refGuided_sc_cluster_guided/data/LP00714.refguided.withGeneSymbols.gtf"}),
}
for s in ["BT474", "K562", "HG002", "UHRR"]:
    DATASETS[s] = dict(kind="bulk", bam=f"{BI}/{s}_E0_merged_other_sorted.bam", genome=BULK_GENOME, baseline=1,
        clusters=None, site_prefix=None,
        gtfs={("denovo", "basic"): f"{REPO}/DENOVO_ID/CellLines/CL_{s}_E0_human/processed_prog_results/{TOOL}.gtf",
              ("refGuided", "basic"): f"{REPO}/REF_Guided/CellLines/CL_{s}_E0_human/processed_prog_results/{TOOL}.gtf"})

def splice_hash(chrom, strand, introns):
    """LRAA's hash (util/gtf_to_isoform_exon_splice_patterns.py), over the 1-based intron list."""
    s = "{}:({}){}".format(chrom, strand, [tuple(x) for x in introns])
    return hashlib.blake2s(s.encode(), digest_size=11).hexdigest()

ATTR = re.compile(r'(\S+) "([^"]*)"')
def attrs(s):
    d = collections.defaultdict(list)
    for k, v in ATTR.findall(s):
        d[k].append(v)
    return d

def read_exons(gtf):
    """transcript_id -> dict(chrom, strand, exons sorted, TSS flag, symbol)"""
    tx, flag = {}, {}
    op = gzip.open if gtf.endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9: continue
            if p[2] == "transcript":
                m = re.search(r'transcript_id "([^"]+)"', p[8]); t = re.search(r'TSS "([^"]+)"', p[8])
                if m and t: flag[m.group(1)] = t.group(1) == "True"
            elif p[2] == "exon":
                m = re.search(r'transcript_id "([^"]+)"', p[8])
                r = tx.setdefault(m.group(1), dict(chrom=p[0], strand=p[6], exons=[]))
                r["exons"].append((int(p[3]), int(p[4])))
    for k, r in tx.items():
        r["exons"].sort()
        r["tss_flag"] = flag.get(k, False)
        r["symbol"] = k.split("^")[0] if "^" in k else ""
    return tx

def intron_chain(exons):
    return tuple((exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1))

def chain_str(c):
    return ",".join(f"{a}-{b}" for a, b in c)
def parse_chain(s):
    return tuple(tuple(int(x) for x in p.split("-")) for p in s.split(","))
PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# Tool-neutral candidate sets (bin/make_pool_gtf.py): for each bulk library, the union of the de novo models of the
# paper's tools and GENCODE, scored by the same read and atlas evidence. Keeps the tool comparison from favouring LRAA.
for s in ["BT474", "K562", "HG002", "UHRR"]:
    DATASETS[s + "_pool"] = dict(DATASETS[s], gtfs={("denovo", "basic"): os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out", "pool", f"{s}.pool.gtf.gz")})
