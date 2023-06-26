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
    
in_dir=main_dir+'data/leipzig/15636.fa/20230331/nii/Part2/'
proc_dir=main_dir+'processing/15636.fa/loraks/'
out_dir=main_dir+'derivatives/15636.fa/loraks/'

folders = ['ernst_kp_mtflash3d_v1q_0p5_sagCS3_LORAKS/', 'pdw_kp_mtflash3d_v1s_0p6CS6_LORAKS/', 't1w_kp_mtflash3d_v1s_0p6CS6_LORAKS/']

for folder in folders:
    mag_files = sorted(glob(in_dir+folder+'*_rec-loraks_*part-mag.nii.gz'))
    phs_files = sorted(glob(in_dir+folder+'*_rec-loraks_*part-phase.nii.gz'))

    den_files = []
    for mag_file in mag_files:
        mag_file = mag_file.replace(in_dir+folder,'').replace('.nii.gz','_den-lcpca110.nii.gz')
        den_files.append(mag_file)

    for phs_file in phs_files:
        phs_file = phs_file.replace(in_dir+folder,'').replace('.nii.gz','_den-lcpca110.nii.gz')
        den_files.append(phs_file)

    denoised = nighres.intensity.lcpca_denoising(image_list=mag_files, phase_list=phs_files, 
                                  ngb_size=4, stdev_cutoff=1.10,
                                  min_dimension=0, max_dimension=-1,
                                  unwrap=True, rescale_phs=False, process_2d=False,
                                  save_data=True, overwrite=True, output_dir=proc_dir+folder,
                                  file_names=den_files)
