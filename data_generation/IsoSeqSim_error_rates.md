# IsoSeqSim error rates

The Arabidopsis and mouse IsoSeqSim simulations used the same four nominal
sequencing-error configurations.

| Dataset suffix | Substitutions (`--es`) | Deletions (`--ed`) | Insertions (`--ei`) | Total error rate | Approximate accuracy |
|---|---:|---:|---:|---:|---:|
| `e0000` | 0.0% | 0.0% | 0.0% | **0.0%** | **100.0%** |
| `e0016` | 0.4% | 0.6% | 0.6% | **1.6%** | **98.4%** |
| `e0050` | 1.7% | 1.1% | 2.2% | **5.0%** | **95.0%** |
| `e0085` | 2.7% | 1.8% | 4.0% | **8.5%** | **91.5%** |

The total error rate is the sum of the configured substitution, deletion, and
insertion rates. The approximate accuracy is `100% - total error rate`. These
are nominal simulation settings; realized rates in the generated reads may
differ slightly because errors are introduced stochastically.

## Sources

- `01. isoseqsim_simulation_Arabidopsis.ipynb`
- `03.isoseqsim_simulation_mouse.ipynb`
