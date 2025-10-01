#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=malann@ccf.org
#SBATCH --job-name=compute_psd
#SBATCH --time=2-02:30:00
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --partition=defq
#SBATCH --output=/mnt/beegfs/malann/codes/Pain_project/logs/%A_%a.out
#SBATCH --error=/mnt/beegfs/malann/codes/Pain_project/logs/%A_%a.err
#SBATCH --array=0 #-33  # Update to 0-23 if using full subject-module list

# Load your environment
#module load python/3.11 # Adjust version as needed
#python -m venv /mnt/beegfs/malann/codes/Pain_project/pain_env

# Activate pre-existing virtual env
source /mnt/beegfs/malann/codes/Pain_project/pain_env/bin/activate
#pip install --upgrade pip
#pip install --no-cache-dir numpy scipy matplotlib mne mne-connectivity nilearn joblib


# List of subjects
SUBJECTS=(
    "PASC001_cb_L_R_good"
    "PASC003_ck_R_L_good"
    "PASC004_tc"
    "PASC005_af_R_L_good"
    "PASC007_ee"
    "PASC009_rg_R_L_good"
    "PASC010_bk_R_L_good"
    "PASC011_ls_R_L_good"
    "PASC012_kk1_R_L_good"
    "PASC013_va_L_R_good"
    "PASC014_sl_L_R_good"
    "PASC015_rj_R_L_good"
    "PASC018_mm_L_R_good"
    "PASC019_bo_L_R_good"
    "PASC020_lj_R_L_good"
    "PASC022_po_R_L_good"
    "PASC023_bj_R_L_good"
    "PASC024_ko_L_R_good"
    "PASC025_hr_L_R_good"
    "PASC027_yj_R_L_good"
    "PASC028_kp_L_R_good"
    "PASC029_aa_L_R_good"
    "PASC030_kk2_L_R_good"
    "PASC032_pm"
    "PASC033_tb_L_R_good"
    "PASC034_mk_R_L_good"
    "PASC035_sk_R_L_good"
    "PASC036_ss"
    "PASC037_gr_L_R_good"
    "PASC038_ge_R_L_good"
    "PASC039_fk_R_L_good"
    "PASC040_tj_L_R_good"
    "PASC041_aa"
)

SUBJECTS=('FASC004_tk_L')

SUBJECT=${SUBJECTS[$SLURM_ARRAY_TASK_ID]}
echo "Running subject: $SUBJECT"

# Run your Python script
python3 /mnt/beegfs/malann/codes/Pain_project/scripts/dSPM_Control.py "$SUBJECT"



