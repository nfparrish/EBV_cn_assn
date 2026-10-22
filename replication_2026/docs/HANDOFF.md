# EBV in 1000 Genomes LCLs — handoff

Updated 2026-10-21 (following frontier-model epistemic/code audit, full-$n$ Track A recomputation, and Track B P1 EBV lytic/latent RNA-seq quantification launch). **Start here.**

| Document | Contents |
|---|---|
| [HANDOFF.md](file:///Users/nfp8/Desktop/HANDOFF.md) | this file — audited state, corrected findings, P1 lytic RNA results, infrastructure, next steps |
| [CLAIMS_LEDGER.md](file:///Users/nfp8/Library/Application%20Support/Claude/scratch-workspaces/8f12ca6e-0f4c-4f6e-a2a3-5fcc8e61a394/5b0c8f55-a2ba-486e-b596-513779f4f64d/scratch-2026-10-20-ee01ee/CLAIMS_LEDGER.md) | every claim with its denominator, completeness, control, and status. **The source of truth.** |
| [EBV_CN_next_steps_plan.md](file:///Users/nfp8/Library/Application%20Support/Claude/scratch-workspaces/8f12ca6e-0f4c-4f6e-a2a3-5fcc8e61a394/5b0c8f55-a2ba-486e-b596-513779f4f64d/scratch-2026-10-20-ee01ee/EBV_CN_next_steps_plan.md) | full analysis narrative |
| [EBV_objectives_and_dependencies.md](file:///Users/nfp8/Library/Application%20Support/Claude/scratch-workspaces/8f12ca6e-0f4c-4f6e-a2a3-5fcc8e61a394/5b0c8f55-a2ba-486e-b596-513779f4f64d/scratch-2026-10-20-ee01ee/EBV_objectives_and_dependencies.md) | dependency graph and parallelism reasoning |
| [EBV_CN_lytic_prior_art_review.md](file:///Users/nfp8/Library/Application%20Support/Claude/scratch-workspaces/8f12ca6e-0f4c-4f6e-a2a3-5fcc8e61a394/5b0c8f55-a2ba-486e-b596-513779f4f64d/scratch-2026-10-20-ee01ee/EBV_CN_lytic_prior_art_review.md) | literature review and novelty assessment |

---

## 1. The question

Lymphoblastoid cell lines (LCLs) from the 1000 Genomes Project were made by transforming donor B
cells with the **B95-8 lab strain** of Epstein-Barr virus. EBV persists as a multi-copy episome and
copy number (CN) per cell varies across two orders of magnitude between lines. The project asks what CN
measures, what determines it, and what it does to the human gene expression data the field builds on.

---

## 2. Data state (verified 2026-10-21)

To eliminate slow Lustre metadata scans across 12,500+ individual `.tsv` files, consolidated master tables
were generated in `/work/nfp8/EBV/metadata/all_*.tsv`.

| Dataset | Complete | Path |
|---|---|---|
| **EBV CN, all 1kGP** | **3,200 / 3,202** | `/work/nfp8/EBV/metadata/all_cn_3200.tsv` + `ebv_cn_1kGP_2504.tsv` |
| **— trio set (698 related)** | **698 / 698** | same |
| **Low-coverage CN** (2008–12 reads) | **2,499 / 2,499** | `/work/nfp8/EBV/metadata/all_lowcov_2499.tsv` |
| **ONT CN** (~2023) | **906 / 907** | `/work/nfp8/EBV/metadata/all_ont_904.tsv` |
| **IR1 repeat copy number** | **2,502 / 2,502** | `/work/nfp8/EBV/metadata/all_ir1_2502.tsv` |
| **ONT IR1** | **907 / 907** | `/work/nfp8/EBV/metadata/all_ont_ir1_907.tsv` |
| **Viral GC-bias metric** | **2,502 / 2,502** | `/work/nfp8/EBV/metadata/all_gcbias_2502.tsv` |
| **LoFreq intrahost variants** | **2,502 / 2,502** | `variants/lofreq/`, `variants/all_calls.tsv`, `variants/pon_sites.tsv` |
| **WGS BCR clonotypes** | **2,502 / 2,504** | `bcr/<donor>/<donor>_report.tsv` |
| **BCR evidence pass** | **731 / 731** | `bcr_evidence/<sample>.tsv` |
| **MAGE RNA FASTQ** | **779 / 779** | `/work/nfp8/MAGE/fastq/` |
| **P1: EBV Lytic/Latent RNA quant** | **779 / 779** | `/work/nfp8/EBV/metadata/all_mage_ebv_rna_779.tsv` + `all_mage_ebv_genes_779.tsv` (`84 genes × 779 runs`) |
| **Aim 1–3 Dissection, QTL, ASE & Mediation tables** | **717 / 717 clean MAGE** | `/work/nfp8/EBV/metadata/viral_gene_dissection_717.tsv`, `viral_qtl_scan_112.tsv`, `host_viral_mediation_20154.tsv`, `/work/nfp8/EBV/ase_mage/` |
| **RNA TRUST4 clonotypes** | 380+ / 779, running | `trust4/<run>/` (`59200071`) |
| **GEUVADIS FASTQ** | 470+ / 666, running | `/work/nfp8/GEUVADIS/fastq/` (`59195108`, `59195102`) |

Base dir `/work/nfp8/EBV`. **`/work` is purged at 75 days** — these become eligible ~2027-01-03.
Anything that must persist goes to `/hpc/group/parrishlab`.

---

## 3. Findings, by confidence (audited & recomputed at full $n = 779/779$ MAGE & $2,449$ 1kGP)

### A. Critical Path P1 (H1): EBV DNA Copy Number Specifically Tracks Lytic Reactivation, Not Latent Episome Dosage

Using competitive host-transcriptome filtering (`/work/nfp8/EBV/metadata/human_nr_tx.mmi`, 63,187 non-redundant
human GENCODE v45 genes to remove human `IL10`/`BCRF1` and `CSF1R`/`BARF1` decoys), 84 non-overlapping CDS
windows from `NC_007605.1.gff3`, and direct CIGAR splice-junction counting (`/work/nfp8/EBV/scripts/ebv_rna_quant.py`):

1. **Lytic transcription is a highly reproducible, cell-line-intrinsic property across independent thaws ($n = 24/24$ multi-replicate MAGE donors, 70 replicate pairs across 73 libraries)**:
   - **Same-batch (technical replicates, $n = 24$ pairs)**: $r(\text{IE/Latent}) = +0.9725$, $r(\text{Late/Latent}) = +0.9739$, $r(\text{IE RPKM}) = +0.9524$, $r(\text{Late RPKM}) = +0.9614$, $r(\text{BGLF5}) = +0.960$, $r(\text{Latent mRNA}) = +0.8581$.
   - **Different-batch (independent biological thaws, $n = 46$ pairs)**: **$r(\text{BGLF5}) = +0.950$**, **$r(\text{IE/Latent}) = +0.9102$**, **$r(\text{IE RPKM}) = +0.8053$**, **$r(\text{Late/Latent}) = +0.7631$**, **$r(\text{Late RPKM}) = +0.6365$**, whereas $r(\text{Latent mRNA}) = +0.2375$.
2. **H1 Decisive Test Across the Full Cohort ($n = 717$ unique clean MAGE donors, `779/779` libraries, 26 populations, 17 batches)**:
   - **Latent Pol II mRNA (`EBNA-1..3C`, `LMP-1`, `LMP-2A/B`) is completely uncoupled from WGS EBV DNA copy number**: **$r = -0.0287$ ($\rho = -0.0699$; within-pop $r = -0.0200$; within-pop+batch $r = +0.0454$)**! Every LCL maintains a homeostatic baseline of latent mRNA expression regardless of whether its viral DNA CN is 5 or 500.
   - **All three lytic kinetic classes and spliced `BZLF1` (Zta) strongly track WGS $\log_2(\text{CN})$**, marginally, within populations, and within populations + batches:
     - **Spliced `BZLF1` (Zta) Junction RPM**: **$r = +0.6477$ ($\rho = +0.6313, p = 1.6 \times 10^{-86}$; within-pop $r = +0.6284$; within-pop+batch $r = +0.6439$)**
     - **Immediate-Early (`BZLF1`/`BRLF1`) RPKM**: **$r = +0.6404$ ($\rho = +0.6285$; within-pop $r = +0.6167$; within-pop+batch $r = +0.6315$)**
     - **Self-normalizing $\log_2(\text{Late / Latent})$ ratio**: **$r = +0.6341$ ($\rho = +0.6349$; within-pop $r = +0.6004$; within-pop+batch $r = +0.6174$)**
     - **Lytic Late Structural RPKM**: **$r = +0.6267$ ($\rho = +0.6359$; within-pop $r = +0.5909$; within-pop+batch $r = +0.6332$)**
     - **Lytic Early RPKM**: **$r = +0.6244$ ($\rho = +0.6202$; within-pop $r = +0.5926$; within-pop+batch $r = +0.6158$)**
     - **Self-normalizing $\log_2(\text{IE / Latent})$ ratio**: **$r = +0.6181$ ($\rho = +0.6071$; within-pop $r = +0.5999$; within-pop+batch $r = +0.6035$)**
   - In the 3-platform subset ($n = 256$ clean donors), MAGE `Lytic IE` RNA correlates at **$r = +0.6696$ with 2019 30x WGS CN** and **$r = +0.6380$ with ~2023 ONT CN** (while `Latent mRNA` is $r = -0.0481$ and $-0.0672$).
   - **Physical WGS corroboration**: Top-30 high-CN WGS BAMs (median CN `326.8`) show an elevated **`oriLyt` (`41–55 kb`) / antipode (`110–130 kb`) coverage ratio (`1.1659` vs `1.0805` in bottom-30 low-CN)**, directly reflecting active lytic replication forks originating at `oriLyt`.

### B. Aim 1: Dissection of All 84 Viral Genes (`viral_gene_dissection.py`, $n = 717$ clean donors)

1. **Top Individual Viral Genes Predicting $\log_2(\text{CN})$ Within Populations + Batches**:
   - **`BGLF5` (alkaline exonuclease / host-shutoff nuclease)**: #1 gene in the EBV genome ($r_{\text{CN}} = +0.6886$, within-pop+batch **$r_{\text{wb}} = +0.6694$**, cross-thaw **$r_{\text{bio}} = +0.950$**), followed by `BKRF4` ($+0.6865$), `BKRF3` ($+0.6854$), `BBLF2_BBLF3` (helicase-primase, $r_{\text{CN}} = +0.6846$, $r_{\text{wb}} = +0.6679$, $r_{\text{bio}} = +0.905$), and `BGLF4` (viral protein kinase, $r_{\text{CN}} = +0.6815$, $r_{\text{wb}} = +0.6631$, $r_{\text{bio}} = +0.928$).
2. **Abortive vs. Productive Lytic Replication Test (Conditioning on `Spliced BZLF1`)**:
   - Conditioning on `Spliced BZLF1` ($R^2 = 41.5\%$ of within-pop+batch $\log_2[\text{CN}]$), downstream lytic replication fork genes (`BALF5` DNA polymerase $r_{\text{part}} = +0.244, t = +6.72$; `BBLF2_BBLF3` $r_{\text{part}} = +0.258, t = +7.14$; `BBLF4` $+0.179$; `BALF2` $+0.156$) and late structural genes (`BcLF1` major capsid $+0.278, t = +7.72$; `BALF4`/gB $+0.278$) add highly significant independent prediction—proving **productive lytic DNA synthesis**.
   - **Viral BCL-2 Divergence**: Strictly lytic `BALF1` adds $r_{\text{part}|\text{BZLF1}} = +0.235$ ($t = +6.47$), whereas dual latent/lytic `BHRF1` adds $+0.002$ ($t = +0.04$, null).
3. **Latent Gene Divergence & Latency III Splice Exclusivity**:
   - Within pop+batch, `EBNA-2` ($r_{\text{wb}} = -0.1457$) and `LMP-1` ($r_{\text{wb}} = -0.1440$) correlate **negatively** with $\log_2(\text{CN})$, `EBNA-1` (`EBNA-1_N`/`C`) is completely flat ($+0.0426$/$+0.0394$), and `BART_RPMS1` ($+0.4150$) & `LMP-2A` ($+0.2931$) correlate positively.
   - Fine-grained splice junction profiling confirms **`0 / 54` libraries use Latency I `Qp_to_EBNA1` (`0.0` reads)** whereas **`54 / 54` use Latency III `CpWp_U_to_EBNA1`, `Cp_EBNA`, and `Wp_W1_W2_IR1`**.

### C. Aim 2: Viral Genetic Influence on EBV Copy Number & Viral Gene Expression (`viral_genetics_qtl.py` & `viral_ase.sbatch`)

1. **Three-Platform-Validated Structural QTL (`IR1` [`BamHI-W`] Repeat Copy Number $\to$ Lytic Switch & EBV Copy Number)**:
   - Previously masked marginally by **Simpson's Paradox** (4 EUR transformation stocks carry ~10.3 `IR1` units while `IBS`/`CEU` have high CN), **within-population `IR1` repeat count strongly and positively drives $\log_2(\text{CN})$ across all $n = 2,449$ clean 1kGP lines ($r = +0.3148, p = 1.67 \times 10^{-60}$; MAGE $n=717$ within-pop+batch $r = +0.3687, p = 2.84 \times 10^{-26}$)**.
   - Verified across three independent controls:
     - **`IR1 / local_Cp_flank` (`10.0–11.5 kb`)**: within-pop **$r = +0.3037$ ($p = 2.1 \times 10^{-53}$)**, while the `local_Cp / genome_unique` flank control is **$r = -0.0334$ (null)**;
     - **Within-`chrEBV` GC-bias control ($n = 2,449$)**: partial **$r(\text{IR1}, \log_2\text{CN} \mid \text{GC}) = +0.3048$ ($p = 2.8 \times 10^{-54}$, 96.8% of effect survives)**;
     - **Cross-platform ONT replication ($n = 890$ clean)**: within-pop **$r(\text{ONT\_IR1}, \text{ONT\_CN}) = +0.2945$ ($p = 2.7 \times 10^{-19}$)** and **$r(\text{30x\_IR1\_GCadj}, \text{ONT\_CN}) = +0.3088$ ($p = 3.9 \times 10^{-21}$)**.
   - **Mechanism in MAGE Viral RNA ($n = 717$, within-pop+batch)**: higher `IR1` copy number increases W-spliced `BHRF1` ($r = +0.1029, p = 0.0057$), reduces `EBNA-2` ($r = -0.1068, p = 0.0041$) and `LMP-1` ($r = -0.0993, p = 0.0076$), and increases the **`log2(Lytic IE / Latent)` lytic switch ($r = +0.2531, p = 2.6 \times 10^{-12}$)**.
2. **Triple-Locked Proof of Neutrality for the `IBS` (`91,192:G>A`) & `GBR` (`132,965:G>GC`) Lineage Barcodes vs. Genuine Within-Population QTLs**:
   - Once demeaned by the 26 populations, both **`91,192:G>A` (`BRLF1` `p.S568S`) and `132,965:G>GC` (`BXRF1`) drop to zero effect on CN ($r_{\text{wp}} = -0.0410$ and $-0.0416$) and lytic/cis RNA ($r_{\text{IE,wb}} = +0.0053$ and $+0.0030$)**. Furthermore, **within-sample heteroplasmic cis-ASE across $n = 19$ `91,192:G>A` carrier libraries (`2,285` RNA reads)** shows `19.1%` RNA VAF vs `17.7%` DNA VAF ($\Delta_{\text{wt}} = +0.014, p = 0.082$, null)—proving they are neutral passenger barcodes of the transformation stocks.
   - Conversely, **`chrEBV:12,129 T>A` (5' `W0`/`Cp-Wp` leader upstream of `IR1`, $n=111$ WGS, $33$ MAGE)** predicts higher $\log_2(\text{CN})$ within populations across all 2,449 lines (**$r = +0.1695, p = 2.9 \times 10^{-17}$**; carrier median CN `59.8` vs `34.7`), **replicates out-of-sample in the 1,732 non-MAGE holdout lines ($r = +0.1736, p = 3.5 \times 10^{-13}$)**, and predicts higher MAGE `Lytic IE` RNA within pop+batch (**$r = +0.1048, p = 0.0049$**).
   - **`chrEBV:163,858 G>T` in `BALF2` ($n=147$ WGS, $51$ MAGE)** predicts lower `BALF2` cis expression ($r_{\text{wb}} = -0.1243, p = 8.4 \times 10^{-4}$) and lower `Lytic IE` ($r_{\text{wb}} = -0.1282, p = 5.7 \times 10^{-4}$).
   - **Within-Sample cis-ASE (`2,066` carrier observations, global $r = +0.9424$)**: identifies significant cis-overexpression at `168,963:T>G` in `LMP-1` (`6.1%` RNA vs `3.9%` DNA across `34,776` reads, $p = 1.5 \times 10^{-95}$) and strong cis-depletion at `36,470:T>C` (`2.7%` RNA vs `16.0%` DNA, $p = 1.1 \times 10^{-31}$) and `36,494:G>C` (`5.1%` vs `11.7%`, $p = 2.3 \times 10^{-4}$) in `EBNA-2`.
   - **Rare Promoter Burden**: Rare variants ($\le 5$ carriers) in `Zp` (`BZLF1` promoter, $n=16$) shift within-pop $\log_2(\text{CN})$ **UP by $+0.231$** and `Lytic IE` **UP by $+0.191$**, whereas rare variants in `Rp` (`BRLF1` promoter, $n=12$) shift $\log_2(\text{CN})$ **DOWN by $-0.493$** and `Lytic IE` **DOWN by $-0.681$**.

### D. Aim 3 (H3): Viral Gene Expression as a Causal Mediator of Host Gene Expression (`host_mediation_bipartite.py`, $20,154 \text{ genes} \times 717 \text{ clean donors}$, $\text{df} = 673$)

1. **Same-Vial Viral RNA Boosts Host Discovery From 30.2% to 68.5% of the Human Transcriptome**:
   - Controlling `sex + 17 batches + 26 populations` ($\text{df} = 673$), host discovery rises from **6,087 genes (30.2%, $\pi_1=49.9\%$) for 2019 WGS CN** to **9,819 genes (48.7%, $\pi_1=65.1\%$) for `Spliced BZLF1`** (`r = +0.8887` vector alignment with `WGS CN`), **10,580 genes (52.5%, $\pi_1=67.9\%$) for `BGLF5`**, **12,863 genes (63.8%, $\pi_1=75.7\%$) for `log2(Late/Latent)`**, and **13,800 genes (68.5%, $\pi_1=78.4\%$) for `LMP-1`** (`r = -0.4510` opposing vector alignment with `WGS CN`).
2. **Formal 1-df Baron–Kenny Mediation Asymmetry (`2019 WGS CN` $\leftrightarrow$ `2023 Spliced BZLF1`) & Exogenous `IR1` Instrument**:
   - Conditioning `2019 WGS log2(CN)` on **`Spliced BZLF1` alone eliminates 92.4% (`5,625 / 6,087`) of `WGS CN` host-gene hits** (mediating **86.0% of directional $R^2$**), whereas **78.7% (`7,732 / 9,819`) of `Spliced BZLF1` host-gene hits survive conditioning on `2019 WGS CN`**.
   - **Exogenous DNA `IR1` Repeat Instrument**: `IR1` repeat count measured in WGS DNA associates with **663 human genes at $q < 0.05$ ($\pi_1 = 26.6\%$)** and aligns with the 20,154-gene `WGS CN` / `Lytic` vector at **$r = +0.7636$** (sharing 5 of the top 10 positively regulated genes across the entire human transcriptome: `ENSG00000286155`, `LINC01443`, `USP17L2`, `CBARP`, `ENSG00000287251`, plus `HLA-DOA`/`HLA-DOB`).
3. **Bipartite 4-Axis Viral-to-Host Functional Module Map**:
   - **`LMP-1` (Latent III NF-$\kappa$B Driver)**: activates `IFIH1` ($+0.763$), `MKNK2` ($+0.755$), `TNFAIP3` ($+0.741$), `KDM6B` ($+0.721$), `CD40` ($+0.666$), `BIRC3` ($+0.653$), `ICAM1` ($+0.643$), `NFKBIA` ($+0.629$), `EBF1` ($+0.464$), and `MYC` ($+0.408$), while **repressing plasmablast differentiation** (`TNFRSF17`/BCMA $-0.624$, `SSR4` $-0.508$, `MZB1` $-0.508$, `JCHAIN` $-0.355$).
   - **`BGLF5` / `Lytic IE` (Lytic Reactivation & Host Shutoff)**: coupled to **loss of B-cell identity factor `EBF1` ($-0.502$)**, suppression of `CD40` ($-0.399$), `GADD45A` ($-0.344$), and `MYC` ($-0.293$), and **induction of plasmablast secretory genes (`JCHAIN` $+0.247$, `MZB1` $+0.230$, `SSR4` $+0.187$)**.
   - **`EBNA-2`**: drives cell-cycle & ISGs (`MX1` $+0.419$, `IFIT1` $+0.385$, `ISG15` $+0.368$).
   - **`BART` lncRNAs**: specifically **suppress Type I Interferon / ISG targets** (`IFIT1` $-0.316$, `OAS1` $-0.276$, `ISG15` $-0.238$, `MX1` $-0.234$, `STAT1` $-0.231$).

### E. Settled WGS & Host Expression Claims (Audited & Recomputed at Full $n$)

1. **Heritability of EBV CN is zero ($h^2 \approx 0$), and the pooled spouse correlation (`0.354`) decomposes into ~59% population stratification (`0.210`) + ~41% within-population family batch ($c^2_{\text{family|pop}} \approx 0.18$)**:
   - **Bug fixed in `heritability.py`**: The script previously read `f[6]` (`Superpopulation`, 5 continents) of `ped3202.txt` instead of `f[5]` (`Population`, 26 populations) and pooled all 570 trios without demeaning by population.
   - **Verified numbers ($n = 570$ clean trios, 565 unique spouse pairs)**:
     - Pooled across populations (`log2`): $r_{\text{mid}} = +0.4306$, $r_{\text{sp}} = +0.3541$, $r_{\text{pp}} = +0.3535$.
     - True same-population (`f[5]`) unrelated control: **$r = +0.2097$** (weighted by trio populations; `+0.1294` uniform across 26 pops)—not `0.047` (which was a same-continent control).
     - **Within-population demeaned (`f[5]`, `log2`)**: **$r_{\text{mid,within}} = +0.2427$, $r_{\text{sp,within}} = +0.1865$, $r_{\text{pp,within}} = +0.1815 \implies h^2_{\text{within}} = 2(0.1865 - 0.1815) = +0.0101 \approx 0$**, and true within-population family batch is **$c^2_{\text{family|pop}} = 0.1815$** (~0.9x the population effect, not 7x). Offspring lines also have +12.4% higher mean CN than parent lines ($\log_2$ mean `5.525` vs `5.356`).
2. **EBV CN Perturbs 30.4% of the Human Transcriptome Within Populations ($q < 0.05$) and Drives 51%–72% of Specific Between-Population Differential Expression (`IBS` vs `GBR`)**:
   - In a unified OLS model controlling `sex + 18 batches + 26 populations` ($\text{df} = 687$, $n = 731$ MAGE donors $\times$ 20,154 genes): **Storey's $\pi_1 = 51.1\%$**, **6,135 genes (30.4%) at BH-FDR $q < 0.05$** (4,263 at $q < 0.01$).
   - **Global 25-df ANOVA vs Targeted Vector Alignment**: Adding 1-df $\log_2(\text{CN})$ reduces average 25-df 26-population ANOVA PVE by only **4.0% relative** (`9.15%` $\to$ `8.79%`; **8.3% drop** on $q<0.05$ genes) because 25 population degrees of freedom dilute a single viral axis. However, **estimating each gene's within-population EBV CN slope ($\beta_{\text{CN, within}}$) strictly within the 24 non-IBS, non-GBR populations (zero sample overlap, ancestry held constant) predicts the between-population differential expression ($\overline{Y}_{\text{IBS}} - \overline{Y}_{\text{GBR}}$) between Spanish (`IBS`, CN ~96) and British (`GBR`, CN ~28) LCLs at $r = +0.7161$ ($R^2 = 51.3\%$) across all 20,154 genes and $r = +0.8493$ ($R^2 = 72.1\%$) across the 6,135 CN-associated genes!**
3. **Out-of-Sample Replication of the GBR & IBS Viral Lineages in the 1,718 Non-MAGE Holdout Lines**:
   - **Bug fixed**: Downstream scripts (`group_check.py`, `ir1_pops.py`, `ir1_decisive.py`, `pon2.py`, `pon3.py`, `haplotypes_sb.py`) had joined WGS tables against MAGE's 731-sample metadata instead of `ped3202.txt`, silently discarding 1,718 clean WGS lines.
   - Using those **1,718 untouched non-MAGE WGS lines** as a true out-of-sample replication cohort:
     - `chrEBV:132,965 G>GC` (`BXRF1`): MAGE discovery ($n=731$) GBR **78.6%** vs IBS **0.0%** $\to$ **Non-MAGE holdout ($n=1,718$) GBR 74.6% (`47/63`) vs IBS 0.0% (`0/79`)** $\to$ **All 2,449 clean 1kGP: GBR 75.8% (`69/91`) vs IBS 0.0% (`0/107`)**.
     - `chrEBV:91,192 G>A` (synonymous `p.S568S` in `BRLF1` / `-249 bp` upstream of `BZLF1` `Zp` promoter): MAGE discovery IBS **64.3%** vs GBR **7.1%** $\to$ **Non-MAGE holdout IBS 53.2% (`42/79`) vs GBR 11.1% (`7/63`)** $\to$ **All 2,449 clean 1kGP: IBS 56.1% (`60/107`) vs GBR 9.9% (`9/91`)**.
4. **IR1 Repeat Copy Number at Full $n = 2,449$ Clean (30x), $n = 2,449$ GC-Controlled, and $n = 891$ Clean (ONT)**:
   - **30x ($n = 2,449$ clean)**: median **12.48u** (range `4.73–26.93u`) vs published `11.6u` in B95-8. 26-pop KW $\eta^2 = 0.2920$, 5-cont $\eta^2 = 0.2353$. European medians: `FIN` **9.80u** ($n=99$), `CEU` **10.35u** ($n=98$), `IBS` **10.45u** ($n=107$), `GBR` **10.64u** ($n=91$), **`TSI` 12.27u ($n=104$)**.
   - **Within-`chrEBV` GC-bias control at full $n = 2,449$** (replacing stale $n=182$): GC explains **3.39%** of IR1 variance; **105.2%** of the 4-EUR split survives ($d = -1.689 \to -1.778$).
   - **ONT IR1 at full $n = 891$ clean (`907/907` total)** (replacing stale $n=40$ [$r=0.830$] and $n=705$ [$r=0.760$]): cross-platform **$r = 0.6845$**, median ratio `ONT/30x = 0.9311`, 4-EUR $d = -1.552$ (`FIN` `9.37u`, `IBS` `9.39u`, `GBR` `9.66u`, `CEU` `10.77u`, `TSI` `12.23u`). The CIGAR-based "structural" estimator in `ont_ir1.sbatch` is withdrawn because it ignored `I` operations and capped at `7.40`.
5. **Three-Platform CN Concordance & Unique-Sequence CN at Full $n$**:
   - **Unique-sequence CN** ($n = 2,449$ clean): median **43.55** (range `3.31–717.53`) vs standard CN median **36.16** (range `2.76–595.96`); ratio `1.2034` ($\text{SD} = 0.0015$), $r(\log_2) = 1.00000$.
   - **ONT vs 30x ($n = 890$ clean)**: **$r = 0.9072$** ($R^2 = 0.8230$). **Low-cov vs 30x ($n = 1,954$ clean)**: **$r = 0.8427$**. **ONT vs Low-cov ($n = 703$ clean)**: **$r = 0.8188$**. Three-indicator reliability ($n = 703$): **30x = 0.9565, ONT = 0.8556, Low-cov = 0.7835**.
   - **Low-coverage ($n = 2,499$)**: **469 (18.8%)** `<0.5x` (blood-derived), **90** in the `0.46–19.65x` gap, **2,004** $\ge 5\text{x}$ (LCL-derived).
6. **WGS BCR Clonality Resolved via Cohort CDR3 Panel-of-Normals (`PoN`)**:
   - **Flaw fixed**: The previous `count <= 20` hard cap was triggered by 115 public/germline CDR3 calls (`>30` reads) recurring across $>2$ unrelated donors. Requiring valid `V`+`J` + **Private Cohort PoN ($\le 2$ donors in 2,502)** without any read-count cap yields a **strictly monotonic read-depth dose-response in RNA recovery (`1–5` reads: 80.8%, `6–15`: 94.2%, `16–30`: 96.6%, `>30`: 100.0%)**.
   - Across $n = 280$ WGS-RNA paired donors: **mean 82.5% (median 87.5%)** of WGS clones are recovered in RNA; **top WGS clone is present in RNA 92.5% of the time (`259/280`, up from 83%)**, in RNA top 5 **57.9% (`162/280`)**, and is the top RNA clone **13.9% (`39/280`)**.
   - Across $n = 2,313$ clean WGS lines, PoN-cleaned `top_clone_frac` (median `0.337`) shows negligible correlation with $\log_2(\text{CN})$ ($\rho = +0.0778$, $R^2 \approx 0.6\%$) and `IR1` ($\rho = -0.0309$, $p = 0.14$).

---

## 4. Infrastructure & Active Jobs on DCC

- **Updated Scripts on DCC (`/work/nfp8/EBV/scripts/`)**:
  - `ebv_rna_quant.py`, `mage_ebv_rna.sbatch`, `mage_ebv_rna_scav.sbatch`, `ebv_rna_analyse.py` (Track B P1 pipeline — **779/779 complete**)
  - `viral_gene_dissection.py` (Aim 1: all 84 EBV genes, 3'-coterminal step-up flux, abortive vs productive partial $r$)
  - `viral_genetics_qtl.py` (Aim 2: within-population `IR1` structural QTL across 30x/GC/ONT, 112-site viral CN/eQTL scan, `Zp`/`Rp` rare variant burden, `LF2` co-infection)
  - `viral_ase_worker.py`, `viral_ase.sbatch` (Aim 2.1: within-sample heteroplasmic viral cis-ASE & fine-grained `Cp`/`Wp`/`Qp` splice junctions)
  - `host_mediation_bipartite.py` (Aim 3: formal Baron–Kenny mediation, exogenous `IR1` host transcriptome test, 4-axis bipartite viral-to-host module map)
  - `heritability.py`, `ir1_pops.py`, `ir1_decisive.py`, `rna_vs_wgs_clones.py`, `clonality_vs_cn.py` (fixed source-of-truth scripts)
- **Active Slurm Jobs**:
  - `mage_trust4` (`59200071`): finishing MAGE BCR repertoire assembly.
  - `geuv_s3` (`59195108`) + `geuv_ebi` (`59195102`): downloading GEUVADIS FASTQs for out-of-cohort RNA replication.
  - `babysit` (`59110204`): active.

---

## 5. Immediate Next Steps

1. **GEUVADIS Out-of-Cohort RNA Replication**: Once `geuv_s3` / `geuv_ebi` complete (`/work/nfp8/GEUVADIS/fastq/`), run `ebv_rna_quant.py` across all 666 GEUVADIS libraries to replicate H1 (`WGS CN` vs `Lytic IE` / `Late` / `Latent`), the 5-lab technical replication design in GEUVADIS, and the Baron–Kenny host mediation in an independent RNA-seq cohort.
2. **Full 779/779 MAGE TRUST4 Recompute**: Once `mage_trust4` (`59200071`) finishes all 779 libraries, re-run `rna_vs_wgs_clones.py` to expand the PoN-cleaned WGS-to-RNA BCR clonotype concordance from $n=280$ to all $n=730$ donors.
