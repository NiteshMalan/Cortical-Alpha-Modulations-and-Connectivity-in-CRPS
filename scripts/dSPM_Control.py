import sys
subject_id = sys.argv[1]  # Take subject name from command line argument

import mne
import gc
from scipy.fft import next_fast_len
import numpy as np

dic = {
    'PASC001':'14741', 'PASC003':'14744', 'PASC005':'14926', 'PASC009':'14829',
    'PASC010':'s10016a', 'PASC011':'14983', 'PASC012':'14938', 'PASC013':'14824',
    'PASC014':'14785', 'PASC015':'14875', 'PASC018':'14922', 'PASC019':'14880',
    'PASC020':'14995', 'PASC022':'s14919', 'PASC023':'14990', 'PASC024':'15005',
    'PASC025':'15028', 'PASC027':'s15076', 'PASC029':'s15256', 'PASC030':'s15162',
    'PASC033':'s15185', 'PASC034':'s15248', 'PASC035':'s15304', 'PASC037':'s15503',
    'PASC038':'s15359', 'PASC039':'s15459', 'PASC040':'s15473'
}

Control_list_no_MRI = ["PASC004_tc", "PASC007_ee", "PASC032_pm", "PASC036_ss", "PASC041_aa", 'PASC028_kp_L_R_good', 'FASC004_tk_L']
Control_avoid = ["PASC014_sl_L_R_good"]

print(f"Subject: {subject_id}", flush=True)

# Handling filenames and paths
if subject_id in ['PASC009_rg_R_L_good']:
    num = '3'
else:
    num = '2'

if subject_id in ["PASC012_kk1_R_L_good", "PASC030_kk2_L_R_good"]:
    val = subject_id[7:11]
else:
    val = subject_id[7:10]

spont = 'SPONT_' + num + val + '_raw_quat_tsss.fif'
name = '/data_clean_SPONT' + num + '-raw.fif'
empty_room_name = '/empty_room_' + num + '-raw.fif'
src_num = 4

# Determine data path and subject MRI directory
if subject_id in Control_list_no_MRI:
    if subject_id in  ['PASC028_kp_L_R_good', 'FASC004_tk_L']:
        data_path = '/mnt/beegfs/malann/codes/Pain_project/data/PASC_MEG_data/' + subject_id
    else:
        data_path = '/mnt/beegfs/malann/codes/Pain_project/data/PASC_MEG_data/Resting state only/' + subject_id

    sub_name = subject_id[0:7]
    subjects_dir = '/mnt/beegfs/malann/codes/Pain_project/data/mne_data/MNE-fsaverage-data'
    subject = "fsaverage"

elif subject_id == "PASC001_cb_L_R_good":
    data_path = '/mnt/beegfs/malann/codes/Pain_project/data/PASC_MEG_data/' + subject_id + '/1702806968/'
    sub_name = subject_id[0:7]
    mri_data_path = "/mnt/beegfs/malann/codes/Pain_project/data/Freesurfer Files/Controls_BEM"
    subjects_dir = mri_data_path + '/subject' + dic[sub_name]
    subject = "sample"

else:
    data_path = '/mnt/beegfs/malann/codes/Pain_project/data/PASC_MEG_data/' + subject_id
    sub_name = subject_id[0:7]
    mri_data_path = "/mnt/beegfs/malann/codes/Pain_project/data/Freesurfer Files/Controls_BEM"
    subjects_dir = mri_data_path + '/subject' + dic[sub_name]
    subject = "sample"

print("Reading raw data...", flush=True)
raw_fname = data_path + '/' + name
raw = mne.io.read_raw_fif(raw_fname, verbose='error')

if subject_id == "PASC014_sl_L_R_good":
    raw.crop(tmin=75, tmax=175)
else: 
    raw.crop(tmin=15, tmax=175)

#raw.resample(1000)
raw.pick_types(meg=True)
sfreq = raw.info['sfreq']

print("Reading empty-room data...", flush=True)

empty_raw_fname = data_path + '/' + empty_room_name
empty_room = mne.io.read_raw_fif(empty_raw_fname,
    allow_maxshield=True
)
empty_room.pick_types(meg=True)
noise_cov = mne.compute_raw_covariance(empty_room)

print("Reading forward solution...", flush=True)
forward_sol = f"/mnt/beegfs/malann/codes/Pain_project/data/PASC_MEG_data/Forward_files/{sub_name}_SPONT_2_ico{src_num}_raw-fwd.fif"
fwd = mne.read_forward_solution(forward_sol)

print("Computing inverse operator...", flush=True)
inverse_operator = mne.minimum_norm.make_inverse_operator(
    raw.info, fwd, noise_cov, loose=0.2, depth=0.8
)
del noise_cov, fwd
gc.collect()

snr = 3
lambda2 = 1.0 / snr**2
n_fft = 2048 * 16 # 16

print("Computing source PSD...", flush=True)
stc_psd = mne.minimum_norm.compute_source_psd(
    raw, inverse_operator, lambda2=lambda2, method='dSPM', overlap=0.5,
    fmin=1, fmax=40.0, n_fft=n_fft,  pick_ori=None, label=None, nave=1, pca=True,prepared=False,
    dB=False, return_sensor=False, verbose=False
)

del inverse_operator, raw
gc.collect()

# Morph to fsaverage
if subject_id not in Control_list_no_MRI:
    print("Morphing to fsaverage...", flush=True)
    sample = "fsaverage"
    fname_fsaverage_src = f"/mnt/beegfs/malann/codes/Pain_project/data/mne_data/MNE-fsaverage-data/{sample}/bem/fsaverage-ico-{src_num}-src.fif"
    src_to = mne.read_source_spaces(fname_fsaverage_src)
    morph = mne.compute_source_morph(
        stc_psd, subject_from="sample", subject_to=sample,
        src_to=src_to, subjects_dir=subjects_dir
    )
    stc_psd = morph.apply(stc_psd)

# Save PSD
#save_path = f"/mnt/beegfs/malann/codes/Pain_project/results/Control_Source_PSDs/{sub_name}_{spont[0:8]}psd_dSPM"
save_path = f'/mnt/beegfs/malann/codes/Pain_project/results/Control_Source_PSDs/{sub_name}_SPONT_2_psd_dSPM'

print(f"Saving PSD to: {save_path}", flush=True)
stc_psd.save(save_path, overwrite=True)
del stc_psd
gc.collect()

print(f"Finished processing: {subject_id}", flush=True)
