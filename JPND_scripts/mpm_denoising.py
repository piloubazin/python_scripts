import nighres
import numpy
import nibabel
from glob import glob
import os
import subprocess
import math
import scipy.ndimage
import json

main_dir='/home/Public/jpnd/'
    
in_dir=main_dir+'data/leipzig/traveling_phantom_7T/mpm_data/'
proc_dir=main_dir+'processing/traveling_phantom_7T/'
out_dir=main_dir+'derivatives/traveling_phantom_7T/'

# input names
prefix = 's2022-12-15_19-21-'

# main processing steps: denoising, T2* fitting 
mpm_files = sorted(glob(in_dir+'*/'+prefix+'*.nii.gz'))

den_files = []
for mpm_file in mpm_files:
    den_file = mpm_file.replace('.nii.gz','_c110.nii.gz')
    den_files.append(den_file)

denoised = nighres.intensity.lcpca_denoising(image_list=mpm_files, phase_list=None, 
                                  ngb_size=4, stdev_cutoff=1.10,
                                  min_dimension=0, max_dimension=-1,
                                  unwrap=False, rescale_phs=False, process_2d=False,
                                  save_data=True, overwrite=True, output_dir=proc_dir,
                                  file_names=den_files)

mpm_files = sorted(glob(in_dir+'*/'+prefix+'*.nii.gz'))
den_files = []
for mpm_file in mpm_files:
    den_file = mpm_file.replace('.nii.gz','_c105.nii.gz')
    den_files.append(den_file)

denoised = nighres.intensity.lcpca_denoising(image_list=mpm_files, phase_list=None, 
                                  ngb_size=4, stdev_cutoff=1.05,
                                  min_dimension=0, max_dimension=-1,
                                  unwrap=False, rescale_phs=False, process_2d=False,
                                  save_data=True, overwrite=True, output_dir=proc_dir,
                                  file_names=den_files)

den_files = []
for mpm_file in mpm_files:
    den_file = mpm_file.replace('.nii.gz','_c100.nii.gz')
    den_files.append(den_file)

denoised = nighres.intensity.lcpca_denoising(image_list=mpm_files, phase_list=None, 
                                  ngb_size=4, stdev_cutoff=1.00,
                                  min_dimension=0, max_dimension=-1,
                                  unwrap=False, rescale_phs=False, process_2d=False,
                                  save_data=True, overwrite=True, output_dir=proc_dir,
                                  file_names=den_files)

