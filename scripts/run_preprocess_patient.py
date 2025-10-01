#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import mne

# --- Get subject name from SLURM argument ---
# Usage: python run_preprocess_patients.py PASP001_md1_NL_AR
if len(sys.argv) < 2:
    raise ValueError("Please provide a patient name as the first argument.")

sub = sys.argv[1]
print(f"Processing patient: {sub}")

# --- Determine subject-specific parameters ---
num = '2'
if sub in ['PASP021_bj1_AL_NR', 'PASP019_bj2_AR_NL']:
    val = sub[7:11]
else:
    val = sub[7:10]

spont = f'SPONT_{num}{val}_raw_quat_tsss.fif'

# --- Data path ---
data_path = f'/mnt/beegfs/malann/codes/Pain_project/data/PASP_MEG_data/{sub}'
raw_fname = os.path.join(data_path, spont)

raw = mne.io.read_raw_fif(raw_fname, preload=True)
info = raw.info

# --- Filtering ---

raw.filter(l_freq=1, h_freq=40, n_jobs=2)
#meg_picks = mne.pick_types(raw.info, meg=True, eeg=True)
#raw.notch_filter(freqs=60, picks=meg_picks)

# Special notch for PASP016_ca_NR_AL
if sub == 'PASP016_ca_NR_AL':
    raw.notch_filter(freqs=13.7, picks=meg_picks)

# --- ICA ---
ica = mne.preprocessing.ICA(n_components=0.95,  method='fastica', max_iter=200, random_state=97)
picks = mne.pick_types(raw.info, meg=True, eeg=False, eog=False,
                       stim=False, exclude='bads')
ica.fit(raw, picks=picks,decim=3, reject=dict(mag=5e-12, grad=5000e-13))

muscle_idx_auto, _ = ica.find_bads_muscle(raw)
ecg_indices, _ = ica.find_bads_ecg(raw, method='ctps', ch_name='ECG001',  threshold="auto") #method="correlation",
try:
    eog_indices, eog_scores = ica.find_bads_eog(raw)
except RuntimeError as e:
    print('Error with EOG detection, msg: {}'.format(e))
    eog_indices, eog_scores = [], None

# Limit to maximum 2 components each
max_ecg = 2
max_eog = 2

if len(ecg_indices) > max_ecg:
    ecg_indices = ecg_indices[:max_ecg]
if len(eog_indices) > max_eog:
    eog_indices = eog_indices[:max_eog]

ica.exclude = ecg_indices + eog_indices + muscle_idx_auto


ica.save(os.path.join(data_path, f'ica{num}-raw.fif'), overwrite=True)

# --- Apply ICA and save cleaned data ---
empty_room = mne.io.read_raw_fif(
    '/mnt/beegfs/malann/codes/Pain_project/data/PASC_MEG_data/empty room 03212025/SPONT_1_re_raw.fif',
    allow_maxshield=True
)
empty_room.load_data()
#empty_room.resample(1000)
empty_room.crop(0,87)
ica.apply(raw)
ica.apply(empty_room)
raw.save(os.path.join(data_path, f'data_clean_SPONT{num}-raw.fif'), overwrite=True)
empty_room.save(os.path.join(data_path, f'empty_room_{num}-raw.fif'), overwrite=True)

print(f"Finished processing {sub}")
