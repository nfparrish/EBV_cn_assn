# Turn-Key Replication Verification Report

Generated automatically by `uv run replication_2026/reproduce_all.py` from `replication_2026/data/`.

| Section | Claim Verified | Recomputed Value | Status |
|---|---|---|---|
| 1. WGS CN | EBV CN range & median (n=2,502 all / n=2,449 clean) | `range=2.76–596.0, med_all=35.69, med_clean=36.16` | **PASS** |
| 1. WGS CN | Within-population trio heritability (n=570 clean trios) | `pooled r_sp=+0.3541/r_pp=+0.3535 | within-pop r_sp=+0.1865/r_pp=+0.1815 -> h²=+0.0101` | **PASS** |
| 1. WGS CN | Three-platform CN concordance (30x vs ONT vs Low-cov) | `r(30x,ONT)=+0.9072 (n=890), r(30x,Low)=+0.8427 (n=1954), r(ONT,Low)=+0.8188 (n=703)` | **PASS** |
| 1. WGS CN | Phase 3 low-coverage blood (<0.5x) vs LCL (>=5x) split | `n_total=2499, blood(<0.5x)=469 (18.8%), LCL(>=5x)=2004` | **PASS** |
| 1. WGS CN | PoN-cleaned WGS BCR clonality vs log2(CN) (n=2,340) | `n=2340, median top_clone_frac=0.432, Spearman rho=+0.0791` | **PASS** |
| 2. H1 RNA | Replicate thaw reproducibility (24 donors, 70 pairs) | `same-batch n=24 r(IE/Lat)=+0.9725 | diff-batch n=46 r(IE/Lat)=+0.9102 vs r(Latent)=+0.2378` | **PASS** |
| 2. H1 RNA | Full-cohort H1 (n=717 clean donors, 779/779 libraries) | `Latent r=-0.0286 (wb -0.0582) | Spliced BZLF1 r=+0.6476 (wb +0.6511) | Lytic IE r=+0.6402 (wb +0.6271)` | **PASS** |
| 3. Aim 1 | Top viral gene BGLF5 & cross-thaw repeatability | `BGLF5 r_CN=+0.6886, r_wb=+0.6694, r_bio=+0.9504` | **PASS** |
| 3. Aim 1 | Productive replication fork & late capsid partial r | BZLF1 | `BALF5(Pol)=+0.2439, BcLF1(Capsid)=+0.2775, BALF1=+0.2354 vs BHRF1=+0.0016` | **PASS** |
| 4. Aim 2 | Three-platform within-population IR1 -> log2(CN) structural QTL | `30x r_wp=+0.3148, GC-controlled r=+0.3048 (96.8%), ONT r_wp=+0.2945 (n=890)` | **PASS** |
| 4. Aim 2 | Out-of-sample replicated 12,129:T>A CN-QTL vs neutral barcodes | `12129:T>A r_wp=+0.1695 (holdout +0.1736) | 91192:G>A r_wp=-0.0410 | 132965:G>GC r_wp=-0.0416` | **PASS** |
| 4. Aim 2 | Within-sample cis-ASE at 91,192:G>A & Latency III vs I splicing | `91192:G>A RNA=541/2633 (0.205) vs DNA=0.192 | Qp_to_EBNA1=0/59 vs CpWp_U_to_EBNA1=59/59` | **PASS** |
| 5. Aim 3 | Host discovery boost & exogenous IR1 instrument alignment | `q<0.05: WGS_CN=6,087, BZLF1=9,819, BGLF5=10,580, LMP1=13,800, IR1_DNA=663 (vector r=+0.7636)` | **PASS** |
| 5. Aim 3 | Zero-sample-overlap vector alignment (24 pops -> IBS vs GBR) | `All 20,154 genes r=+0.7086 (R²=50.2%) | 6,087 CN-sig genes r=+0.8444 (R²=71.3%)` | **PASS** |
| 6. Aim 4 | Viral-conditioned trio heritability (n=570 clean trios) | `M1 r_pp=+0.1815 -> M4 r_pp=+0.1119 (-38.3%) | M4 h²=+0.0578 (95% CI [-0.122, +0.228])` | **PASS** |
| 6. Aim 4 | Genome-wide 7.62M host GWAS & DNMT3B lytic switch candidate | `GWAS P<5e-8: PopAdj=1 -> ViralCond=0 (max PVE=1.05%) | DNMT3B IE/Lat P=1.63e-06 (q=0.0220)` | **PASS** |
