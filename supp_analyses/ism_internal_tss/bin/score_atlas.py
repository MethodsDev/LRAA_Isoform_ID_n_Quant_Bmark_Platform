"""The CAGE/refTSS distance helper, shared."""
import numpy as np, pandas as pd
from common import *
class Atlas:
    def __init__(self, bed):
        d = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2, 3, 4, 5], names=["c", "s", "e", "id", "score", "st"], dtype={"c": str})
        self.g = {k: (v.s.values, v.e.values, v.id.values, v.score.values) for k, v in
                  ((k, v.sort_values("s")) for k, v in d.groupby(["c", "st"]))}
    def dist(self, c, st, p):
        p = np.atleast_1d(p)
        if (c, st) not in self.g:
            return np.full(len(p), 10**9), np.array([""] * len(p), dtype=object), np.zeros(len(p))
        s, e, ids, sc = self.g[(c, st)]
        i = np.searchsorted(s, p - 1, side="right")
        best = np.full(len(p), 10**9); idx = np.zeros(len(p), dtype=int)
        for off in (-3, -2, -1, 0):
            j = np.clip(i + off, 0, len(s) - 1)
            d = np.where(p - 1 < s[j], s[j] - (p - 1), np.where(p - 1 >= e[j], (p - 1) - e[j] + 1, 0))
            u = d < best
            best = np.where(u, d, best); idx = np.where(u, j, idx)
        return best, ids[idx], sc[idx]
cage, reftss = Atlas(CAGE), Atlas(REFTSS)
