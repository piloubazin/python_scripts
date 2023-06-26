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
    
in_dir=main_dir+'data/leipzig/38714.7a/20230404/nii/'
proc_dir=main_dir+'processing/38714.7a/20230404/nii/'
out_dir=main_dir+'derivatives/38714.7a/20230404/nii/'

folders = ['t1w_kp_mtflash3d_v1s_0p8CS6_LORAKSlp1/',\
           'mtw_kp_mtflash3d_v1s_0p8CS6_LORAKSlp1/',\
           'pdw_kp_mtflash3d_v1s_0p8CS6_LORAKSlp1/']
           
mag_files = []
for folder in folders:
    new_files = sorted(glob(in_dir+folder+'*_rec-loraks_*part-mag.nii.gz'))
    for new_file in new_files:
        mag_files.append(new_file)

den_files = []
for mag_file in mag_files:
    mag_file = mag_file.replace(in_dir,'').replace('.nii.gz','_den-lcpca-m120-n3-joint.nii.gz')
    den_files.append(mag_file)

print('denoised files: '+str(den_files))

denoised = nighres.intensity.lcpca_denoising(image_list=mag_files, phase_list=None, 
                                  ngb_size=3, stdev_cutoff=1.20, use_rmt=False,
                                  min_dimension=0, max_dimension=-1,
                                  unwrap=True, rescale_phs=False, process_2d=False,
                                  save_data=True, overwrite=True, output_dir=proc_dir,
                                  file_names=den_files)
