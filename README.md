# Cortical Alpha Modulations and Connectivity in CRPS

## 🧠 Overview

This repository contains analysis pipelines for investigating **Peak Alpha Frequency (PAF)**, **Alpha Power**, and **functional connectivity** in **Complex Regional Pain Syndrome (CRPS)** using EEG/MEG data.

The workflow integrates:

* MNE source-space analysis
* FreeSurfer head modeling
* Brainstorm conversion
* SPM statistical modeling
* Python-based connectivity analysis
* HPC batch processing

---

## 📂 Repository Structure

```
.
├── jobscripts/
│   ├── run_preprocess_control.sh
│   ├── run_preprocess_patient.sh
│   ├── run_psd_control.sh
│   └── run_psd_patient.sh
│
├── SPM-MATLAB-Codes/
│   ├── read_mne_PAF.m
│   ├── read_mne_PAF_Power.m
│   ├── BS_gifti_PAF_script.m
│   └── BS_gifti_PAF_Power_script.m
│
├── notebooks/
│   ├── Compare_source_PSDs.ipynb
│   ├── Connectivity_Patient.ipynb
│   ├── Connectivity_control.ipynb
│   ├── Forward_solution_Control.ipynb
│   ├── Forward_solution_Patient.ipynb
│   └── freesurfer_head_model.ipynb
│
├── scripts/
│   ├── dSPM_Control.py
│   ├── dSPM_Patient.py
│   ├── run_preprocess_control.py
│   └── run_preprocess_patient.py
│
└── README.md
```

---

## 📁 Folder Description

### `jobscripts/`

SLURM scripts for HPC processing:

* `run_preprocess_control.sh` — preprocessing controls
* `run_preprocess_patient.sh` — preprocessing patients
* `run_psd_control.sh` — PSD computation (controls)
* `run_psd_patient.sh` — PSD computation (patients)

---

### `SPM-MATLAB-Codes/`

#### PAF Analysis

* `read_mne_PAF.m` — Reads MNE `.stc` PAF files into Brainstorm
* `BS_gifti_PAF_script.m` — Converts to SPM and performs PAF–pain correlation

#### Alpha Power Analysis

* `read_mne_PAF_Power.m` — Reads MNE `.stc` power files into Brainstorm
* `BS_gifti_PAF_Power_script.m` — Converts to SPM and performs Power–pain correlation

PAF and Alpha Power pipelines are handled **independently**.

---

### `notebooks/`

* `freesurfer_head_model.ipynb` — FreeSurfer head model creation
* `Forward_solution_Control.ipynb` — Forward solution for controls
* `Forward_solution_Patient.ipynb` — Forward solution for patients
* `Compare_source_PSDs.ipynb` — Source PSD comparison
* `Connectivity_control.ipynb` — Connectivity analysis (controls)
* `Connectivity_Patient.ipynb` — Connectivity analysis (patients)

---

### `scripts/`

Python processing scripts:

* `run_preprocess_control.py` — preprocessing pipeline for controls
* `run_preprocess_patient.py` — preprocessing pipeline for patients
* `dSPM_Control.py` — source reconstruction (dSPM) for controls
* `dSPM_Patient.py` — source reconstruction (dSPM) for patients

---

## 🔬 Analysis Workflow

### 1. Head Model & Forward Solution

1. Create FreeSurfer head model
2. Compute forward solutions for controls and patients

### 2. Preprocessing

1. Run Python preprocessing scripts
2. Execute SLURM job wrappers

### 3. Source Reconstruction

1. Compute dSPM source estimates
2. Export `.stc` files

### 4. Spectral Analysis

1. Compute source PSD
2. Extract alpha band
3. Compute PAF and Alpha Power

### 5. SPM Statistical Analysis

1. Import `.stc` into Brainstorm
2. Convert to GIFTI
3. Run correlation with pain ratings in SPM

### 6. Connectivity Analysis

1. ROI extraction
2. Connectivity matrix computation
3. Group comparison
4. Visualization

---

## 🚀 Running on HPC

Preprocessing:

```bash
sbatch jobscripts/run_preprocess_control.sh
sbatch jobscripts/run_preprocess_patient.sh
```

PSD computation:

```bash
sbatch jobscripts/run_psd_control.sh
sbatch jobscripts/run_psd_patient.sh
```

---

## 📊 Outputs

* FreeSurfer head models
* Forward solutions
* dSPM source estimates
* Source PSD comparisons
* PAF maps
* Alpha Power maps
* Pain correlation statistics
* Connectivity matrices
* Publication-ready figures

---

## 🎯 Features

✔ FreeSurfer head modeling
✔ dSPM source reconstruction
✔ Separate PAF and Power pipelines
✔ SPM pain correlation analysis
✔ Connectivity analysis
✔ HPC-ready SLURM scripts
✔ Publication-quality visualization

---

## 👨‍🔬 Author

Nitesh Malan
