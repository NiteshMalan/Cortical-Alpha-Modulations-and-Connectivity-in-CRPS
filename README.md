Cortical-Alpha-Modulations-and-Connectivity-in-CRPS# Cortical Alpha Modulations and Connectivity in CRPS

## 🧠 Overview

This repository contains analysis pipelines for investigating **cortical alpha-band modulations** and **functional connectivity** in **Complex Regional Pain Syndrome (CRPS)** using EEG/MEG data.

The project includes:

* Alpha-band power analysis
* Peak Alpha Frequency (PAF) estimation
* ROI-based connectivity analysis
* Bayesian statistical testing
* Age/sex matched comparisons
* Publication-quality connectivity visualizations

---

## 📂 Repository Structure

```id="f3wkjv"
.
├── jobscripts/            # SLURM job submission scripts for HPC processing
├── SPM-MATLAB-Codes/      # MATLAB/SPM scripts for preprocessing and source analysis
├── notebooks/             # Jupyter notebooks for exploration and figure generation
├── scripts/               # Python analysis pipelines and utilities
└── README.md
```

### Folder Description

#### `jobscripts/`

Contains batch scripts for running large-scale analyses on HPC clusters:

* Connectivity computation
* Spectral analysis
* Permutation testing
* Parallel subject-level processing

#### `SPM-MATLAB-Codes/`

MATLAB scripts using SPM for:

* Source reconstruction
* ROI extraction
* MEG/EEG preprocessing
* Time-frequency analysis

#### `notebooks/`

Interactive notebooks for:

* Data exploration
* Statistical analysis
* Connectivity visualization
* Figure generation for publications

#### `scripts/`

Core Python pipeline including:

* Connectivity computation
* Bayes factor analysis
* ROI filtering
* Matched cohort comparison
* Heatmap generation

---

## ⚙️ Requirements

Python ≥ 3.9

Main dependencies:

```id="9skpvz"
numpy
scipy
pandas
matplotlib
seaborn
mne
pingouin
scikit-learn
```

Install:

```bash id="r8nfjv"
pip install -r requirements.txt
```

---

## 🔬 Analysis Workflow

1. Preprocessing (SPM / MATLAB)
2. ROI extraction
3. Connectivity matrix computation
4. Statistical comparison (CRPS vs Controls)
5. Bayes Factor calculation
6. Age/Sex matched analysis
7. Visualization and figure generation

---

## 📊 Connectivity Analysis

* ROI-based connectivity matrices
* Upper triangular extraction
* Group comparison
* Bayes Factor (BF10) computation
* Directional effect labeling
* Matched cohort comparison

---

## 🧪 Statistical Interpretation

| BF10    | Evidence          |
| ------- | ----------------- |
| < 1/3   | Evidence for null |
| 1/3 – 3 | Anecdotal         |
| 3 – 10  | Moderate          |
| 10 – 30 | Strong            |
| > 30    | Very strong       |

---

## 🚀 Example Usage

Run Python connectivity analysis:

```bash id="07b82r"
python scripts/run_connectivity_analysis.py
```

Run on HPC:

```bash id="e0xdr3"
sbatch jobscripts/connectivity_job.sh
```

---

## 📈 Outputs

The pipeline generates:

* Connectivity matrices
* Bayes factor tables
* Significant ROI lists
* Publication-ready heatmaps
* Matched vs full cohort comparison figures

---

## 🎯 Features

✔ Alpha-band connectivity analysis
✔ CRPS vs Control comparison
✔ Bayesian inference
✔ Age/sex matched cohort analysis
✔ HPC-ready job scripts
✔ SPM-based source analysis
✔ Publication-ready visualization

---

## 👨‍🔬 Author

Nitesh Malan

---

## 📄 License

MIT License (update if needed)

---

## 🤝 Acknowledgments

* MNE-Python
* SPM
* Pingouin
* Seaborn
