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
#SBATCH --array=0-32  # Fixed indexing to match full subject list

# Activate virtual environment
source /mnt/beegfs/malann/codes/Pain_project/pain_env/bin/activate

# List of subjects
SUBJECTS=('PASP001_md1_NL_AR' 'PASP002_lr_AR_NL' 'PASP003_kb_NR_AL' 'PASP004_fm_NR_AL' 'PASP005_je_AL_NR' 
          'PASP006_ed' 'PASP007_kj_AR_NL' 'PASP008_gv_NR_AL' 'PASP009_kk_AR_NL' 'PASP010_ka' 'PASP011_gc_NL_AR'
          'PASP012_ct_AR_NL' 'PASP013_yd_NL_AR' 'PASP014_bs_NL_AR' 'PASP015_da_NR_AL' 'PASP016_ca_NR_AL'
          'PASP017_md2_AL_NR' 'PASP018_ht' 'PASP019_bj2_AR_NL' 'PASP020_rp_NL_AR' 'PASP021_bj1_AL_NR'
          'PASP022_sa_NL_AR' 'PASP023_ie_NL_AR' 'PASP025_wl_AR_NL' 'PASP026_mj_AL_NR' 'PASP027_me_AL_NR'
          'PASP028_lt_NL_AR' 'PASP029_wr_AL_NR' 'PASP030_hj_NR_AL' 'PASP031_pp_AL_NR' 'PASP032_nd_AL_NR'
          'PASP033_sa2_NL_AR')

SUBJECT=${SUBJECTS[$SLURM_ARRAY_TASK_ID]}
echo "Running subject: $SUBJECT"

python3 /mnt/beegfs/malann/codes/Pain_project/scripts/dSPM_Patient.py "$SUBJECT"
