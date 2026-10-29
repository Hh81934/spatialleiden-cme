# spatialleiden-cme

Stability analysis of SpatialLeiden clustering on MERFISH mouse hypothalamus data using Cross Method Entropy (CME).

This repository contains the analysis notebook and scripts from a Bachelor's internship project. SpatialLeiden was run across 625 parameter combinations, and Cross Method Entropy (CME) was used to map which regions of the tissue are clustered consistently and which are sensitive to parameter choice.

---

## Background

Spatial clustering of spatial transcriptomics data is an exploratory task with no established guidelines. Results depend on several choices:

- **Layer ratio:** how much weight is given to spatial location compared to gene expression
- **Latent resolution:** how finely cells are split by gene expression
- **Spatial resolution:** how finely cells are split by spatial location
- **Random seed:** the random initialization of the algorithm

CME computes an entropy value for every cell across many clustering runs, independent of the cluster labels:

- **Low entropy (green):** the cell is assigned consistently across runs (stable region)
- **High entropy (red):** the cell's assignment changes between runs (unstable region)

---

## Dataset

- **Technology:** MERFISH (Multiplexed Error-Robust Fluorescence In Situ Hybridization)
- **Tissue:** mouse hypothalamus preoptic region
- **Size:** 5,488 cells, 155 genes
- **File:** `MERFISH_0.04.h5ad` (included in this repository)
- **Reference:** Moffitt et al. 2018

---

## Repository Structure

```
spatialleiden-cme/
├── spatialleiden_cme.ipynb      main analysis notebook
├── MERFISH_0.04.h5ad            MERFISH dataset
├── merfish_env.yml              conda environment
├── figures/                     figures shown in the notebook
├── cme_script_original.R        original R CME implementation
├── cme_script_modified.R        R script used for the K=4 vs K=6 control experiment
├── cme_script_fix3bugs.R        corrected R CME implementation
├── cme_script_original.py       Python CME implementation
└── README.md
```

The notebook calls the R and Python scripts by file name, so keep them in the same folder as `spatialleiden_cme.ipynb`.

---

## Getting Started

Clone the repository:

```bash
git clone https://github.com/Hh81934/spatialleiden-cme.git
cd spatialleiden-cme
```

Create and activate the conda environment:

```bash
conda env create -f merfish_env.yml
conda activate merfish_env
```

The environment provides Python 3.11 with scanpy, squidpy, spatialleiden, scipy, pandas, numpy and matplotlib.

The R scripts additionally need R 4.x with the packages `clue` and `optparse`:

```r
install.packages(c("clue", "optparse"))
```

Then open `spatialleiden_cme.ipynb` in JupyterLab or JupyterHub and run the cells in order.

---

## Notebook Overview

1. **Introduction:** background, dataset and objectives
2. **SpatialLeiden vs Leiden:** reproduction of the SpatialLeiden tutorial
3. **Parameter scan:** layer ratio, latent resolution and spatial resolution scans, compared to the anatomical atlas
4. **Nested for loop:** 625 SpatialLeiden runs over all parameter combinations
5. **CME pipeline:** aggregation of runs, selection of base clusterings, entropy maps from the R and Python scripts
6. **Control experiments:** K=4 vs K=6, and label swapping on K=8
7. **Debugging the R script:** investigation of the differences between the R and Python implementations
8. **Biological interpretation:** comparison of CME maps to known hypothalamus structures
9. **Summary and conclusions**

Figures are saved to `figures/` and shown in the notebook through markdown image links, which keeps the notebook file small.

### Parameter ranges used for the 625 runs

| Parameter | Values |
|---|---|
| Layer ratio | 1.4 – 1.8 in 0.1 steps |
| Latent resolution | 0.115 – 0.207 in 0.023 steps |
| Spatial resolution | 0.9 – 1.5 in 0.15 steps |
| Random seed | 0, 42, 123, 456, 789 |

5 × 5 × 5 × 5 = 625 runs, giving cluster counts from K=4 to K=13.

---

## Running the CME Scripts

CME follows a three-step workflow:

1. **Aggregate** all clustering labels into one TSV file (`cme_input_all.tsv`): rows are cells, columns are clustering runs.
2. **Select base clusterings** with a `BC_rankings.tsv` file: the column headers are cluster counts and the values are run names. In this project the first 8 runs per cluster count were used.
3. **Run the CME script** for one cluster count at a time.

The notebook creates `cme_input_all.tsv` and `BC_rankings.tsv` for you.

**Corrected R script:**

```bash
Rscript cme_script_fix3bugs.R \
    -i cme_input_all.tsv \
    -o cme_output_k8.tsv \
    --BC_ranking BC_rankings.tsv \
    --n_bcs 8 \
    -n 8
```

**Python script:**

```bash
python cme_script_original.py \
    -i cme_input_all.tsv \
    -o cme_output_k8.tsv \
    --BC_ranking BC_rankings.tsv \
    --n_bcs 8 \
    -n 8
```

**Options**

| Option | Description |
|---|---|
| `-i` | input TSV with aggregated labels |
| `-o` | output TSV with one entropy value per cell |
| `--BC_ranking` | TSV listing the base clusterings per cluster count |
| `--n_bcs` | number of base clusterings to use |
| `-n` | cluster count K to compute CME for |
| `--no_split` | (Python only) do not split runs by cluster count |

K=4 has only 2 runs in the parameter scan, so it was run without `BC_rankings.tsv` and with `--n_bcs 2`.

---

## Summary of Findings

- CME maps were computed for K=4 to K=13 from 625 SpatialLeiden runs.
- Clustering stability generally decreases as the cluster count increases.
- The original R script and the Python implementation gave different CME maps on identical input. Control experiments (K=4 vs K=6, label swapping) located the difference in the label alignment step of the R script.
- The corrected R script (`cme_script_fix3bugs.R`) matches the Python output closely, with correlation 1.0 for 9 of 10 cluster counts and 0.9997 for K=9.
- Stable regions in the corrected CME maps correspond to well-defined hypothalamus structures, and unstable regions lie mostly at the boundaries between structures.

The implementation differences found in the R script are still being investigated and should be validated further.

---

## References

- Sun et al. — Cross Method Entropy (CME). *[add full citation and DOI]*
- Müller-Bötticher et al. (2025). SpatialLeiden: spatially aware Leiden clustering. *Genome Biology*.
- Moffitt et al. (2018). Molecular, spatial, and functional single-cell profiling of the hypothalamic preoptic region. *Science*. https://doi.org/10.1126/science.aau5324
- Anatomical atlas adapted from Li & Zhou et al. 2022.

---

## Author

**Hosna Haydari**
Bachelor's intern, Cancer Bioinformatics Group, BIH Berlin, 2026

Supervision: Prof. Naveed Ishaque
Discussions and feedback: Alexander Malt
