# EBV Copy Number, Viral Genetics, and Host Transcriptome Mediation in 1000 Genomes / MAGE LCLs (`replication_2026`)

Self-contained turn-key replication package, consolidated master tables, source-of-truth claims ledger, and publication figure suite for dissecting Epstein-Barr virus (EBV) episome copy number (`CN`), viral genetic determinants, lytic vs. latent viral transcription, and causal mediation of the human transcriptome across **3,200 1000 Genomes WGS lines**, **907 Oxford Nanopore (`ONT`) long-read genomes**, **2,499 Phase 3 low-coverage genomes**, and **all 779 / 779 MAGE RNA-seq libraries (`730` WGS-matched donors = `717` clean B95-8 donors + `13` wild-type co-infected donors; `70` technical and independent-thaw replicate pairs)**.

---

## 1. Turn-Key Replication Quick Start (<15 Seconds Locally)

All consolidated master tables (`data/*.tsv`) are bundled directly in this repository so that **every statistical test, table, and multi-panel figure can be reproduced deterministically on any laptop or workstation without requiring access to raw BAM/FASTQ files or the Duke Compute Cluster (DCC)**:

```bash
# Run the complete verification suite (14/14 ledger assertions) + regenerate all 6 multi-panel figures:
uv run replication_2026/reproduce_all.py

# Or regenerate only the 6 publication figures (PNG 300 DPI + vector PDF):
uv run replication_2026/scripts/generate_figures.py
```

Running `reproduce_all.py` writes the machine-verified output report to [`REPLICATION_VERIFICATION_REPORT.md`](REPLICATION_VERIFICATION_REPORT.md) and saves all 6 multi-panel figures to [`figures/`](figures/).

To re-run the upstream BAM/CRAM/FASTQ extraction and Slurm arrays on Duke DCC (`/work/nfp8/EBV` / `/hpc/group/parrishlab/EBV_1kGP_MAGE`):

```bash
# Re-run all statistical models from consolidated tables on DCC:
bash replication_2026/scripts/run_dcc_pipeline.sh --from-master-tables

# Submit full Slurm job arrays from raw 1kGP CRAMs, ONT BAMs, and MAGE FASTQs:
bash replication_2026/scripts/run_dcc_pipeline.sh --submit-slurm
```

---

## 2. Publication Figure Suite (`figures/`)

| Figure | File (`.png` & `.pdf`) | Summary |
|---|---|---|
| **Figure 1** | [`figures/fig1_wgs_architecture_and_heritability.png`](figures/fig1_wgs_architecture_and_heritability.png) | **Multi-Platform EBV Copy Number Architecture, Zero Host Heritability ($h^2 \approx 0$), & BCR Clonality Independence**: (A) 26-population $\log_2(\text{CN})$ distribution ($n=2,449$ clean 1kGP); (B) 30x Illumina vs. ~2023 ONT long-read concordance ($r=+0.9072, n=890$) and 3-platform reliability triangle; (C) Within-population trio heritability decomposition ($n=570$ clean trios: $r_{\text{sp,within}}=+0.1865$ vs. spouse $r_{\text{pp,within}}=+0.1815 \implies h^2=+0.010 \approx 0$); (D) PoN-cleaned WGS BCR clonality dose-response in RNA (`80.8%` $\to$ `100.0%`) and independence from $\log_2(\text{CN})$ ($\rho=+0.0778, R^2=0.6\%$). |
| **Figure 2** | [`figures/fig2_lytic_vs_latent_rna_h1.png`](figures/fig2_lytic_vs_latent_rna_h1.png) | **Decisive H1 Test — EBV DNA Copy Number Tracks Productive Lytic Reactivation, Not Latent Episome Dosage ($n=779$ MAGE Libraries)**: (A) Cross-thaw biological repeatability ($n=24$ donors, $70$ pairs: `BGLF5` $r_{\text{bio}}=+0.950$, `log2(IE/Latent)` $r_{\text{bio}}=+0.910$ vs. `Latent mRNA` $r_{\text{bio}}=+0.238$); (B) Joint scatter of 2019 WGS $\log_2(\text{CN})$ vs. `Latent Pol II mRNA` ($r=-0.0287$, null) and `Lytic IE RPKM` ($r=+0.6404$); (C) Invariance across marginal, within-population, and within-population+batch models (`Spliced BZLF1` $r=+0.6477$); (D) 100% Latency III (`Cp`/`Wp`) vs. `0 / 59` Latency I (`Qp_to_EBNA1`) splice isoform architecture. |
| **Figure 3** | [`figures/fig3_viral_gene_dissection_84.png`](figures/fig3_viral_gene_dissection_84.png) | **Aim 1 — Dissection of All 84 EBV Genes & Productive Lytic Replication**: (A) Genome-wide ranking of all 84 non-overlapping EBV gene windows against within-pop+batch $\log_2(\text{CN})$, led by host-shutoff nuclease `BGLF5` ($r=+0.6886, r_{\text{wb}}=+0.6694$); (B) Cross-thaw stability ($r_{\text{bio}}$) vs. CN coupling; (C) Abortive vs. productive partial correlation conditioning on `Spliced BZLF1` (`BALF5` Pol $+0.244$, `BBLF2/3` $+0.258$, `BcLF1` capsid $+0.278$, strictly lytic `BALF1` $+0.235$ vs. dual `BHRF1` $+0.002$); (D) 3'-coterminal poly(A) step-up initiation flux. |
| **Figure 4** | [`figures/fig4_viral_genetics_ir1_qtl_and_ase.png`](figures/fig4_viral_genetics_ir1_qtl_and_ase.png) | **Aim 2 — Viral Genetic Determinants (`IR1` Structural QTL, `12,129:T>A` CN-QTL, Lineage Neutrality, & Within-Sample cis-ASE)**: (A) Three-platform-validated within-population `IR1` (`BamHI-W`) repeat copy number QTL ($r=+0.3148, p=1.7\times 10^{-60}$ in 30x; $r=+0.3048$ GC-controlled; $r=+0.2945$ in ONT; $r=+0.3687$ in MAGE); (B) Mechanism of `IR1` repeat number in MAGE RNA (`BHRF1` $+0.103$, lytic switch $+0.253$, `EBNA-2` $-0.107$); (C) 112-site within-population viral QTL scan contrasting neutral `IBS`/`GBR` barcodes (`91,192:G>A`, `132,965:G>GC`) against out-of-sample-replicated `12,129:T>A` ($r=+0.170$, holdout $r=+0.174$), `163,858:G>T` (`BALF2`), and `Zp`/`Rp` rare variant burden; (D) Within-sample heteroplasmic viral cis-ASE ($r=+0.942$ across 2,100+ SNV observations). |
| **Figure 5** | [`figures/fig5_host_mediation_and_bipartite_map.png`](figures/fig5_host_mediation_and_bipartite_map.png) | **Aim 3 — Causal Mediation of the Human Transcriptome & 4-Axis Bipartite Viral-to-Host Module Map ($20,154\text{ Genes} \times 717\text{ Clean Donors}$)**: (A) Same-vial viral RNA boosts host discovery from `6,087` genes (`30.2%`, `2019 WGS CN`) to `9,819` (`48.7%`, `Spliced BZLF1`), `10,580` (`52.5%`, `BGLF5`), and `13,800` (`68.5%`, `LMP-1`); (B) 1-df Baron–Kenny mediation asymmetry (`Spliced BZLF1` eliminates `92.4%` of `WGS CN` hits and mediates `86.0%` of directional $R^2$, whereas `78.7%` of `Spliced BZLF1` hits survive conditioning on `WGS CN`); (C) Exogenous DNA `IR1` repeat instrument ($663$ human genes at $q<0.05$, aligning with `WGS CN` at $r=+0.7636$); (D) 4-axis bipartite viral-to-host heatmap (`LMP-1`, `BGLF5`, `Spliced BZLF1`, `EBNA-2`, `BART`). |
| **Figure 6** | [`figures/fig6_population_vector_alignment_ibs_gbr.png`](figures/fig6_population_vector_alignment_ibs_gbr.png) | **Out-of-Sample Prediction of Between-Population Differential Expression (`24 Non-IBS/GBR Pops` $\to$ `IBS vs. GBR`)**: (A) Within-population viral slope $\hat{\beta}_{\text{CN}}$ trained strictly in the 24 non-IBS/GBR populations (zero sample overlap) predicts observed `IBS vs. GBR` differential expression across all 20,154 human genes ($r=+0.7086$ to $+0.7161$, $R^2=50.2\%\text{–}51.3\%$) and CN-significant genes ($r=+0.8444$ to $+0.8493$, $R^2=71.3\%\text{–}72.1\%$); (B) Multi-parameter viral model comparison (`R² = 73.8%`) vs. global 25-df ANOVA dilution. |

---

## 3. Repository Structure

- [`reproduce_all.py`](reproduce_all.py) — Master turn-key verification and figure generation script (`uv run replication_2026/reproduce_all.py`).
- [`REPLICATION_VERIFICATION_REPORT.md`](REPLICATION_VERIFICATION_REPORT.md) — Automated verification report generated by `reproduce_all.py`.
- [`docs/CLAIMS_LEDGER.md`](docs/CLAIMS_LEDGER.md) — Audited source-of-truth ledger of every claim, denominator, control, and withdrawn/corrected prior result.
- [`docs/HANDOFF.md`](docs/HANDOFF.md) — Executive synthesis and methods documentation.
- [`data/`](data/) — Consolidated master tables:
  - `all_cn_3200.tsv`, `ebv_cn_1kGP_2504.tsv`, `ped3202.txt` — 30x WGS EBV copy number & pedigree metadata ($n=3,200$).
  - `all_ont_904.tsv`, `all_lowcov_2499.tsv` — Oxford Nanopore ($n=906$) and Phase 3 low-coverage ($n=2,499$) EBV CN.
  - `all_ir1_2502.tsv`, `all_ont_ir1_907.tsv`, `all_gcbias_2502.tsv` — `IR1` (`BamHI-W`) repeat copy number & within-`chrEBV` GC bias.
  - `all_mage_ebv_rna_779.tsv`, `all_mage_ebv_genes_779.tsv` — Competitive decoy-filtered, splice-aware EBV RNA quantification across all 779 MAGE libraries and 84 viral genes.
  - `all_mage_viral_ase_59.tsv`, `all_mage_viral_splices_59.tsv` — Within-sample heteroplasmic viral cis-ASE and fine-grained `Cp`/`Wp`/`Qp` splice junction counts.
  - `viral_gene_dissection_84.tsv` — Aim 1 summary table across all 84 EBV genes.
  - `viral_qtl_scan_112.tsv` — Aim 2 within-population viral CN-QTL and eQTL scan across 112 polymorphic sites (with 1,732-line non-MAGE holdout replication).
  - `host_viral_mediation_20154.tsv` — Aim 3 transcriptome-wide associations, 1-df Baron–Kenny mediation statistics, exogenous `IR1` effects, and out-of-sample `24-pop -> IBS vs GBR` vector alignment across all 20,154 human genes.
  - `bcr_wgs_summary.tsv` — Private Cohort Panel-of-Normals (`PoN`)-cleaned WGS `IGH` clonality summary ($n=2,313$ clean lines).
  - `dennis_ko_lcls_annotated_480.tsv` — Annotated Coriell LCL lookup catalog with multi-platform CN, `IR1`, MAGE lytic/latent RNA, and key viral variants.
- [`scripts/`](scripts/) — Complete collection of Python analysis scripts and Slurm `.sbatch` jobs, including [`generate_figures.py`](scripts/generate_figures.py) and [`run_dcc_pipeline.sh`](scripts/run_dcc_pipeline.sh).
