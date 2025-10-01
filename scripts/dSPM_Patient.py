import sys
subject_id = sys.argv[1]  # Take subject name from command line argument

import mne
import gc
from scipy.fft import next_fast_len
import numpy as np

dic = {
'PASP001':'s15363','PASP002':'s15508','PASP004':'15763','PASP014':'17105',
'PASP017':'15913','PASP019':'16850','PASP020':'17053','PASP021':'16803',
'PASP022':'16190','PASP023':'16408','PASP028':'16591','PASP029':'16547',
'PASP030':'16527','PASP031':'16342','PASP033':'16883'
}


Patient_list_with_MRI = ['PASP001_md1_NL_AR', 'PASP002_lr_AR_NL', 'PASP004_fm_NR_AL', 'PASP014_bs_NL_AR', 'PASP017_md2_AL_NR',
                            'PASP022_sa_NL_AR', 'PASP023_ie_NL_AR', 'PASP028_lt_NL_AR', 'PASP029_wr_AL_NR', 'PASP030_hj_NR_AL',
                            'PASP031_pp_AL_NR','PASP019_bj2_AR_NL' , 'PASP021_bj1_AL_NR', 'PASP020_rp_NL_AR', 'PASP033_sa2_NL_AR']

Patient_avoid = [ 'PASP009_kk_AR_NL','PASP011_gc_NL_AR','PASP019_bj2_AR_NL']
Patient_list = ['PASP005_je_AL_NR', 'PASP015_da_NR_AL']

print(f"Subject: {subject_id}", flush=True)


# Handling filenames and paths
src_num = 4
num = '2'
val = subject_id[7:10]
data_path = f'/mnt/beegfs/malann/codes/Pain_project/data/PASP_MEG_data/{subject_id}'
spont = 'SPONT_' + num + val + '_raw_quat_tsss.fif'
name = '/data_clean_SPONT' + num + '-raw.fif'
empty_room_name = '/empty_room_' + num + '-raw.fif'

src_num = 4


# Determine subject-specific configurations
if subject_id in Patient_list_with_MRI:
    sub_name = subject_id[0:7]
    mri_data_path = "/mnt/beegfs/malann/codes/Pain_project/data/Freesurfer Files/Patients_BEM"
    subjects_dir = f"{mri_data_path}/subject{dic[sub_name]}"
    subject = "sample"
else:
    sub_name = subject_id[0:7]
    subjects_dir = "/mnt/beegfs/malann/codes/Pain_project/data/mne_data/MNE-fsaverage-data"
    subject = "fsaverage"



print("Reading raw data...", flush=True)
raw_fname = data_path + '/' + name
raw = mne.io.read_raw_fif(raw_fname, verbose='error')
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
forward_sol = f"/mnt/beegfs/malann/codes/Pain_project/data/PASP_MEG_data/Forward_files/{sub_name}_SPONT_2_ico{src_num}_raw-fwd.fif"
fwd = mne.read_forward_solution(forward_sol)

print("Computing inverse operator...", flush=True)
inverse_operator = mne.minimum_norm.make_inverse_operator(
    raw.info, fwd, noise_cov, loose=0.2, depth=0.8
)
del noise_cov, fwd
gc.collect()

snr = 3
lambda2 = 1.0 / snr**2
n_fft = 2048 * 16

print("Computing source PSD...", flush=True)
stc_psd = mne.minimum_norm.compute_source_psd(
    raw, inverse_operator, lambda2=lambda2, method='dSPM', overlap=0.5,
    fmin=1, fmax=40.0, n_fft=n_fft,  pick_ori=None, label=None, nave=1, pca=True,prepared=False,
    dB=False, return_sensor=False, verbose=False
)



del inverse_operator, raw
gc.collect()


# Morph to fsaverage
if subject_id in Patient_list_with_MRI:
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
#save_path = f"/mnt/beegfs/malann/codes/Pain_project/results/Patient_Source_PSDs/{sub_name}_{spont[0:8]}psd_dSPM"
save_path = f'/mnt/beegfs/malann/codes/Pain_project/results//Patient_Source_PSDs/{sub_name}_{spont[0:8]}psd_dSPM'

print(f"Saving PSD to: {save_path}", flush=True)
stc_psd.save(save_path, overwrite=True)
del stc_psd
gc.collect()

print(f"Finished processing: {subject_id}", flush=True)
