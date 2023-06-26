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
    
in_dir=main_dir+'data/leipzig/traveling_phantom_7T/ernst_sag/'
proc_dir=main_dir+'processing/traveling_phantom_7T/'
out_dir=main_dir+'derivatives/traveling_phantom_7T/'

# input names
ernst_prefix = 's2022-12-15_19-21-201255-00001-'

# main processing steps: denoising, T2* fitting 
ernst_files = sorted(glob(in_dir+ernst_prefix+'*.nii.gz'))
ernst_jsons = sorted(glob(in_dir+ernst_prefix+'*.json'))

ernst_tes = []
for ernst_json in ernst_jsons:
    #print("reading: "+ernst_json)
    data = json.load(open(ernst_json, mode='r'))
    te = data['acqpar'][0]['EchoTime']
    print("TE = "+str(te))
    ernst_tes.append(te)

qr2 = nighres.intensity.flash_t2s_fitting(image_list=ernst_files,
                                te_list=ernst_tes,
                                save_data=True,overwrite=False,
                                output_dir=proc_dir)              

den_files = []
for ernst_file in ernst_files:
    den_file = ernst_file.replace('.nii.gz','_c100.nii.gz').replace(in_dir,'')
    den_files.append(den_file)

denoised = nighres.intensity.lcpca_denoising(image_list=ernst_files, phase_list=None, 
                                  ngb_size=4, stdev_cutoff=1.00,
                                  min_dimension=0, max_dimension=-1,
                                  unwrap=False, rescale_phs=False, process_2d=False,
                                  save_data=True, overwrite=False, output_dir=proc_dir,
                                  file_names=den_files)

qr2 = nighres.intensity.flash_t2s_fitting(image_list=denoised['denoised'],
                                te_list=ernst_tes,
                                save_data=True,overwrite=False,
                                output_dir=proc_dir)              

comb = nighres.intensity.t2s_optimal_combination(image_list=denoised['denoised'],
                                te_list=ernst_tes,
                                save_data=True,overwrite=False,
                                output_dir=proc_dir)
