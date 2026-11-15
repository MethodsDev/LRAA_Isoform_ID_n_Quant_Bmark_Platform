"""Consequence of starting at S's TSS for L's CDS. cds: genomic [(a, b)] of L's CDS (no stop)."""
import numpy as np
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def tx_pos(exons, strand):
    x = np.concatenate([np.arange(a, b + 1) for a, b in sorted(exons)])
    return x if strand == "+" else x[::-1]
def tx_seq(fa, chrom, exons, strand):
    s = "".join(fa.fetch(chrom, a - 1, b).upper() for a, b in sorted(exons))
    return s if strand == "+" else s.translate(COMP)[::-1]
def s_exons(chain, tss, end3, strand):
    c = sorted(chain)
    if strand == "+":
        ex = [(tss, c[0][0] - 1)] + [(c[i][1] + 1, c[i + 1][0] - 1) for i in range(len(c) - 1)] + [(c[-1][1] + 1, end3)]
    else:
        ex = [(end3, c[0][0] - 1)] + [(c[i][1] + 1, c[i + 1][0] - 1) for i in range(len(c) - 1)] + [(c[-1][1] + 1, tss)]
    return [(a, b) for a, b in ex if b >= a]
def kozak(seq, i):
    """seq index i = A of ATG"""
    if i < 3 or i + 3 >= len(seq): return "", "NA"
    ctx = seq[max(0, i - 6):i + 4]; r = seq[i - 3] in "AG"; g = seq[i + 3] == "G"
    return ctx, "strong" if r and g else ("adequate" if r or g else "weak")
def orf_len(seq, i):
    for j in range(i, len(seq) - 2, 3):
        if seq[j:j + 3] in ("TAA", "TAG", "TGA"):
            return (j - i) // 3, j
    return (len(seq) - i) // 3, None      # no stop in S
EMPTY = dict(protein_class="", aa_total=np.nan, aa_lost=np.nan, aa_remaining=np.nan, kozak=np.nan, kozak_ctx=np.nan, first_ATG_in_frame=np.nan,
             first_ATG_orf_aa=np.nan, first_ATG_kozak=np.nan, new_atg_gpos=np.nan, s_ends_before_stop=np.nan, start_to_tss_nt=np.nan)
def classify(fa, chrom, strand, chain, tss, end3, cds):
    ex = s_exons(chain, tss, end3, strand)
    seq, gp = tx_seq(fa, chrom, ex, strand), tx_pos(ex, strand)
    cp = tx_pos(cds, strand)
    out = dict(EMPTY); n_aa = len(cp) // 3; out["aa_total"] = n_aa
    insq = np.isin(cp, gp)
    if insq[0]:
        out["protein_class"] = "start_kept"
        out["start_to_tss_nt"] = int(np.where(gp == cp[0])[0][0])    # 5'UTR length in S
        return out
    if not insq.any():
        out["protein_class"] = "no_CDS_in_S"; return out
    k = int(np.argmax(insq))
    s0 = int(np.where(gp == cp[k])[0][0])
    out["s_ends_before_stop"] = not insq[-1]
    fr0 = s0 - k                         # S index of CDS base 0 (may be negative)
    c = (k + 2) // 3; hit = None
    while 3 * c < len(cp):
        i = fr0 + 3 * c
        if i + 3 > len(seq): break
        if seq[i:i + 3] == "ATG": hit = c; break
        c += 1
    fa_i = seq.find("ATG")
    if fa_i < 0:
        out["protein_class"] = "no_ATG_in_S"; return out
    forf, fstop = orf_len(seq, fa_i)
    out["first_ATG_orf_aa"] = forf; out["first_ATG_kozak"] = kozak(seq, fa_i)[1]
    in_frame = (fa_i - fr0) % 3 == 0
    out["first_ATG_in_frame"] = in_frame
    if hit is not None:
        i = fr0 + 3 * hit
        ctx, kz = kozak(seq, i)
        out.update(aa_lost=hit, aa_remaining=n_aa - hit, kozak=kz, kozak_ctx=ctx, new_atg_gpos=int(cp[3 * hit]))
        if fa_i == i: out["protein_class"] = "N_truncated"
        elif in_frame and fstop is not None and fstop > i: out["protein_class"] = "N_altered_in_frame"   # in-frame ATG upstream, read through into CDS (novel first exon)
        elif fstop is not None and fstop < i and forf < 30: out["protein_class"] = "N_truncated_after_uORF"
        else: out["protein_class"] = "frame_changed"
    else:
        out["protein_class"] = "frame_changed"
    return out
