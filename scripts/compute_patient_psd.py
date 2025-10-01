import sys
import gc
import mne
import numpy as np
from scipy.fft import next_fast_len

def process_subject(sub):
    try:
        print(f"Processing subject: {sub}")

        src_num = 4
        num = '2'
        val = sub[7:10]
        data_path = f'/mnt/beegfs/malann/PASP_MEG_data/{sub}'
        name = f'/data_clean_SPONT{num}-raw.fif'
        spont = 'SPONT_'+num+val+'_raw_quat_tsss.fif'

        # Define patient list with MRI
        Patient_list_with_MRI = [
            'PASP001_md1_NL_AR', 'PASP002_lr_AR_NL', 'PASP004_fm_NR_AL', 'PASP014_bs_NL_AR', 'PASP017_md2_AL_NR',
            'PASP022_sa_NL_AR', 'PASP023_ie_NL_AR', 'PASP028_lt_NL_AR', 'PASP029_wr_AL_NR', 'PASP030_hj_NR_AL',
            'PASP031_pp_AL_NR','PASP019_bj2_AR_NL' , 'PASP021_bj1_AL_NR', 'PASP020_rp_NL_AR', 'PASP033_sa2_NL_AR'
        ]

        # Check if MRI is available
        if sub in Patient_list_with_MRI:
            subject = "sample"
            subjects_dir = f"/mnt/beegfs/malann/Freesurfer Files/Patients_BEM/subject{dic[sub[:7]]}"
        else:
            subject = "fsaverage"
            subjects_dir = "/mnt/beegfs/malann/mne_data/MNE-fsaverage-data"

        # Load raw MEG data
        raw_fname = f"{data_path}/{name}"
        raw = mne.io.read_raw_fif(raw_fname, preload=True)
        raw.crop(tmin=15, tmax=175)       
        raw.pick_types(meg=True)
        sfreq = raw.info['sfreq']
                        
        # Compute noise covariance from empty room data
        empty_room = mne.io.read_raw_fif('/mnt/beegfs/malann/PASC_MEG_data/empty room 03212025/data_clean_SPONT2-raw.fif', allow_maxshield=True)
        empty_room.pick_types(meg=True)
        noise_cov = mne.compute_raw_covariance(empty_room)

        # Load forward model
        fwd = mne.read_forward_solution(f"/mnt/beegfs/malann/PASP_MEG_data/Forward_files/{sub[:7]}_SPONT_2_ico{src_num}_raw-fwd.fif")

        # Compute inverse operator
        inverse_operator = mne.minimum_norm.make_inverse_operator(
            raw.info, fwd, noise_cov, loose=0.2, depth=0.8
        )

        # Free memory
        del noise_cov, fwd
        gc.collect()

        # Compute source power spectral density (PSD)
        snr = 3
        lambda2 = 1.0 / snr**2
        method = "dSPM"
        n_fft = 2048 * 16  # Adjust for memory limits

        stc_psd = mne.minimum_norm.compute_source_psd(
            raw,
            inverse_operator,
            lambda2=lambda2,
            fmin=1,
            fmax=40,
            n_fft=n_fft,
            dB=False,
            return_sensor=False,
            verbose=False,
        )

        # Free memory
        del inverse_operator, raw
        gc.collect()

        # Morph to fsaverage if MRI available
        if sub in Patient_list_with_MRI:
            src_to = mne.read_source_spaces("/mnt/beegfs/malann/mne_data/MNE-fsaverage-data/fsaverage/bem/fsaverage-ico-4-src.fif")
            morph = mne.compute_source_morph(stc_psd, subject_from="sample", subject_to="fsaverage",
                                             src_to=src_to, subjects_dir=subjects_dir)
            stc_psd = morph.apply(stc_psd)

        # Save results
        stc_psd.save(f'/mnt/beegfs/malann/Patient_Source_PSDs/{sub[:7]}_{spont[:8]}_psd_dSPM', overwrite=True)

        print(f"✅ Finished processing {sub}")

    except Exception as e:
        print(f"❌ Error processing {sub}: {e}")

# Run script with command-line subject argument
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python process_meg.py <subject_name>")
        sys.exit(1)

    subject_name = sys.argv[1]
    process_subject(subject_name)
