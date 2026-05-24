**Single-cell and bulk BCR repertoire analysis pipeline for longitudinal HIV data**

We developed a reproducible bioinformatics pipeline for processing, integrating, and analyzing B-cell receptor (BCR) sequencing data used in the study:

> **Autologous Antibodies Targeting the Silent Face of HIV Envelope Exert Potent Virologic Control**  
> Chun et al., *Nature Medicine* (2026, Submitted)

---

## Overview

The pipeline seeks to investigate autologous antibody responses targeting the HIV envelope silent face harmonizing **bulk BCR**, **10x Genomics single-cell**, and **sorted/cultured B cell** datasets within the Immcantation framework.

---

## Data Availability

The processed datasets used in this study are publicly available at:

**Zenodo:** [https://zenodo.org/records/20045799](https://zenodo.org/records/20045799)

---

# Methods

- **Gene Annotation**: We identified V and J gene segments using **IgBLAST** (v1.22.0) against the IMGT germline database
- **Clonal Clustering**:
  - Sequences were grouped by shared V gene, J gene, and junction length
  - Optimal clustering threshold was determined using the `findThreshold` function from **Shazam** v1.3.1 on length-normalized Hamming distances
  - Clonal assignment was performed with single-linkage hierarchical clustering via `hierarchicalClones` from **SCOPer** v1.4.0
  - Clones were further refined by splitting those with mismatched light chain (LC) V/J genes and resolving LC pairing using **Dowser** v2.4.1
- **Phylogenetic Analysis** (C02 and G11 lineages):
  - Trees were inferred using **Dowser** v2.4.1 and **IgPhyML** v2.0.0 with separate HC/LC partitions (HLP19 model)
  - Evolutionary rate (mutations/codon/day) was calculated as the slope of genetic divergence vs. sample time
  - Significance was assessed via date randomization test (10,000 permutations) in **Dowser** v2.4.1

---

## Repository Structure

```bash
hiv-bcr-analysis/
├── scripts/
│   ├── 01_data_preprocessing/          # Data preparation and integration
│   ├── 02_clonal_analysis/             # Clonal clustering
│   ├── 03_phylogenetic_analysis/       # Tree building & mutation analysis
│   └── pipeline/                       # Main SLURM execution scripts
├── results/                            # Final results (trees, plots, tables)
└── log/                                # Pipeline logs
```

## Pipeline Stages
1. IgBLAST Annotation (`01_run_igblast.sh`)
- Annotates sorted and cultured B-cell sequences using `IgBLAST`

2. Dataset Merging & Preprocessing (`02_run_preprocessing.sh`)
- Integrates bulk BCR, 10x single-cell, and sorted/cultured data into standardized AIRR format

3. Clonal Clustering (`03_run_clonal_clustering.sh`)
- Determines optimal clonal clustering threshold and performs hierarchical clonal assignment with light chain refinement

4. Phylogenetic & Mutational Analysis (`04_run_phylogenetics.sh`)
- Builds lineage trees for C02/G11 lineages, calculates evolutionary rates, and runs date randomization tests.

---

## Usage

```bash
cd /path/to/hiv-bcr-analysis

sbatch scripts/pipeline/01_run_igblast.sh
sbatch scripts/pipeline/02_run_preprocessing.sh
sbatch scripts/pipeline/03_run_clonal_clustering.sh
sbatch scripts/pipeline/04_run_phylogenetics.sh
```
---

## Citation

We’re really glad you found this pipeline helpful!

If you use it in your research, we’d greatly appreciate it if you cite the associated publication:

```bibtex
@article{Chunetal2026,
  author = {Chun, Tae-Wook and Shang, Xiaoran and Galkin, Andrey and Zhang, Xiaozhen and Wu, Sherry and de Assis, Felipe Lopes and Jang, Junseok and Kardava, Lela and Buckner, Clarisa and Atkinson, Ben and Burgos, Eduar Pinzon and Wang, Wei and Melnyk, Mattie and Justement, Jesse and Sewack, Adeline and Shi, Victoria and Baxter, Rebecca and Pozharski, Edwin and Gittens, Kathleen and Sneller, Michael and Blazkova, Jana and Hoehn, Kenneth and Heredia, Alonso and Li, Yuxing and Moir, Susan},
  title = {Autologous Antibodies Targeting the Silent Face of HIV Envelope Exert Potent Virologic Control},
  journal = {Nature Medicine},
  year = {2026},
  note = {Submitted}
}
```

Feel free to reach out if you have any questions or need help adapting the pipeline!