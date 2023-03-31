# This script performs parcellation on the fMRI scans and creates a connectivity
# matrix for each subject to be used for further analysis.

from nilearn import image as nimg
from nilearn import maskers
from nilearn.connectome import ConnectivityMeasure
from nilearn.interfaces.fmriprep import load_confounds
import os
import numpy as np
import pandas as pd
import sys
from sys import argv
from pathlib import Path
from pathlib import Path

print("Reading arguments")
subj=argv[1]
parc=argv[2]
outdir=argv[3]

# example run: python parcellation_trial_oct_11.py 'sub-10209' '/home/HCP_MNI_152.nii' '/home/confound_parc_trial'

#Load separated parcellation
parcellation = nimg.load_img(parc)

#Load TR filelist and obtain TR
os.chdir('/home/')
RT_file = pd.read_table("RT_list.txt", sep = ",")
RT_file = np.array(RT_file)
for i in range(len(RT_file)):
        if RT_file[i,0] == subj: 
            TR = RT_file[i,1]

TR = float(TR)

#Load masks
print("Loading parcellation")
masker = maskers.NiftiLabelsMasker(labels_img=parcellation, standardize=True, 
memory='nilearn_cache', verbose=1, detrend=True, low_pass = 0.1, 
high_pass = 0.01, t_r=TR)

# Set directories
print("Setting up directories")
dir_var = '/home/fmriprep/'
if subj == 'sub-10207' or  subj == 'sub-10209' or subj == 'sub-10211' or subj == 'sub-10212' or subj == 'sub-10213' or subj == 'sub-10215' or subj == 'sub-10255' or subj == 'sub-10807':
    file_f = subj+'_task-rest_run-1_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz'
    file_m = subj+'_task-rest_run-1_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz'
elif subj == 'sub-10159':
    file_f = subj+'_task-rest_run-2_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz'
    file_m = subj+'_task-rest_run-2_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz'
else:
    file_f = subj+'_task-rest_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz'
    file_m = subj+'_task-rest_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz'
func_file = dir_var+subj+'/func/'+file_f
mask_file = dir_var+subj+'/func/'+file_m

#Load functional image
print("Loading functional image")
func_img = nimg.load_img(func_file)

#Apply cleaning, parcellation and extraction to functional data
cleaned_and_averaged_time_series = masker.fit_transform(func_img) 
cleaned_and_averaged_time_series.shape
save_time_series = pd.DataFrame(cleaned_and_averaged_time_series) 

Path(outdir+'/'+subj).mkdir(parents=True, exist_ok=True)
os.chdir(outdir+'/'+subj)
print("Saving timeseries")
save_time_series.to_csv('time_series_'+subj+'.txt', sep=',')

# Check if the number of masker labels equals the number of atlas labels that we have
print("Number of labels", len(masker.labels_))

# Get the label numbers from the atlas
atlas_labels = np.unique(parcellation.get_fdata().astype(int))

# Get number of labels that we have
NUM_LABELS = len(atlas_labels)
print(NUM_LABELS)

correlation_measure = ConnectivityMeasure(kind='correlation')
full_correlation_matrix = correlation_measure.fit_transform([cleaned_and_averaged_time_series])
full_correlation_matrix.shape

save_matrix = full_correlation_matrix[0, :, :]
save_matrix = pd.DataFrame(save_matrix) 
#print(save_matrix)

os.chdir(outdir+'/'+subj)
print("Saving connectivity matrices")
save_matrix.to_csv('cor_mat_'+subj+'.txt', sep=',')

# Extract time-series and matrices with confound regression
print("Rerunning with confound regression")
confounds_simple, sample_mask = load_confounds(
    func_file,
    strategy=["high_pass", "motion", "wm_csf"],
    motion="basic", wm_csf="basic")

print("The shape of the confounds matrix is:", confounds_simple.shape)
print(confounds_simple.columns)

time_series = masker.fit_transform(func_file,
                                   confounds=confounds_simple,
                                   sample_mask=sample_mask)
time_series.shape
save_time_series = pd.DataFrame(time_series)

os.chdir(outdir+'/'+subj)
print("Saving timeseries")
save_time_series.to_csv('Deconfounded_time_series_'+subj+'.txt', sep=',')

correlation_matrix = correlation_measure.fit_transform([time_series])[0]
np.fill_diagonal(correlation_matrix, 0)
save_matrix = pd.DataFrame(correlation_matrix)

os.chdir(outdir+'/'+subj)
print("Saving connectivity matrices")
save_matrix.to_csv('Deconfounded_cor_mat_'+subj+'.txt', sep=',')

# Extract time-series and matrices with confound regression
print("Rerunning with confound regression + GSR")
confounds_simple, sample_mask = load_confounds(
    func_file,
    strategy=["high_pass", "motion", "wm_csf","global_signal"],
    motion="basic", wm_csf="basic",global_signal="basic")

print("The shape of the confounds matrix is:", confounds_simple.shape)
print(confounds_simple.columns)

time_series = masker.fit_transform(func_file,
                                   confounds=confounds_simple,
                                   sample_mask=sample_mask)
time_series.shape
save_time_series = pd.DataFrame(time_series)

os.chdir(outdir+'/'+subj)
print("Saving timeseries")
save_time_series.to_csv('DeconfoundedGSR_time_series_'+subj+'.txt', sep=',')

correlation_matrix = correlation_measure.fit_transform([time_series])[0]
np.fill_diagonal(correlation_matrix, 0)
save_matrix = pd.DataFrame(correlation_matrix)

os.chdir(outdir+'/'+subj)
print("Saving connectivity matrices")
save_matrix.to_csv('DeconfoundedGSR_cor_mat_'+subj+'.txt', sep=',')
