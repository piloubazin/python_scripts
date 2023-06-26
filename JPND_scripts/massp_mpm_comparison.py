import nighres
import numpy
import nibabel
from glob import glob
import os
import subprocess
import math
import scipy.ndimage

main_dir='/home/Public/jpnd/'
    
proc0p6_dir=main_dir+'processing/15636.fa/mpm0p6/'
proc0p7_dir=main_dir+'processing/15636.fa/mpm0p7/'
proc_dir=main_dir+'processing/15636.fa/mpm_comparison/'

# MPM input names
qr1_0p6_file = proc0p6_dir+'mpm0p6_r1_brain.nii.gz'
qr2s_0p6_file = proc0p6_dir+'mpm0p6_r2s_brain.nii.gz'
qpd_0p6_file = proc0p6_dir+'mpm0p6_pd_brain.nii.gz'
massp_0p6_file = proc0p6_dir+'mpm0p6_massp-label.nii.gz'

qr1_0p7_file = proc0p7_dir+'mpm0p7_r1_brain.nii.gz'
qr2s_0p7_file = proc0p7_dir+'mpm0p7_r2s_brain.nii.gz'
qpd_0p7_file = proc0p7_dir+'mpm0p7_pd_brain.nii.gz'
massp_0p7_file = proc0p7_dir+'mpm0p7_massp-label.nii.gz'


if (not os.path.exists(proc_dir)):
    os.makedirs(proc_dir)

print("\n1. MPM coreg\n")

ants = nighres.registration.embedded_antspy_multi(
                            source_images=[qr1_0p6_file,qr2s_0p6_file,qpd_0p6_file],
                            target_images=[qr1_0p7_file,qr2s_0p7_file,qpd_0p7_file],
                            run_rigid=True, run_affine=False, run_syn=False,
                            rigid_iterations=10000,
                            cost_function='MutualInformation', 
                            interpolation='Linear',
                            smooth_mask=0.1,
                            ignore_affine=True, 
                            save_data=True, file_name='mpm0p6-to-mpm0p7',
                            output_dir=proc_dir)

# map labels in both spaces
massp0p6_coreg = nighres.registration.apply_coordinate_mappings(massp_0p6_file,
                            mapping1=ants['mapping'],
                            save_data=True, output_dir=proc_dir)['result']

massp0p7_coreg = nighres.registration.apply_coordinate_mappings(massp_0p7_file,
                            mapping1=ants['inverse'],
                            save_data=True, output_dir=proc_dir)['result']

nighres.statistics.segmentation_statistics(massp_0p6_file, intensity=None, template=massp0p7_coreg,
                            statistics=['Volumes','Volume_difference','Dice_overlap','Dilated_Dice_overlap'], 
                            output_csv='mpm_comparison.csv', output_dir=proc_dir)

nighres.statistics.segmentation_statistics(massp_0p7_file, intensity=None, template=massp0p6_coreg,
                            statistics=['Volumes','Volume_difference','Dice_overlap','Dilated_Dice_overlap'], 
                            output_csv='mpm_comparison.csv', output_dir=proc_dir)