import nibabel as nb
import nighres as nr
import glob
import numpy as np
import os
from nighres.io import load_volume, save_volume


# List of subject IDs to be processed, assuming they all come from the Subcortex database (BIDS version)
#subjects = sorted(glob.glob('./r1maps/*.nii.gz'))
#mni_brain = '/usr/share/fsl/data/standard/MNI152_T1_0.5mm.nii.gz'            
#mni_brain = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/MNI152_T1_0.5mm_brain_d4mm.nii.gz'            
mni_t1 = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/mni_icbm152_t1_tal_nlin_asym_09b_hires_brain_d4mm.nii.gz'
mni_t2 = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/mni_icbm152_t2_tal_nlin_asym_09b_hires_brain_d4mm.nii.gz'
mni_pd = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/mni_icbm152_pd_tal_nlin_asym_09b_hires_brain_d4mm.nii.gz'

# Overall approach: 1) use nifti header + rigid registration to align whole brain and slab data, 
# 2) register whole brain to MNI (rigidly), 3) transfer slab data to MNI space. Coregistration 
# usually works better from whole brain to slab, so we do that and take the inverse mapping. 
# Additional labels for the whole brain data can also be mapped back onto the slab that way.

subjects = [0,1,2,4,8,12,24,31,33,40]

for subject in subjects: 
    subject = str(subject).zfill(3)
    brain_r1 = os.path.join('/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/r1maps/',
                            'sub-'+subject+'_ses-1_acq-wb_mod-r1map_orient-std_brain.nii.gz')
    brain_r2s = os.path.join('/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/r2smap/',
                            'sub-'+subject+'_ses-1_acq-wb_mod-r2starmap_orient-std_brain.nii.gz')
    brain_qsm = os.path.join('/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/qsm/',
                            'sub-'+subject+'_ses-1_acq-wb_mod-qsm_orient-std_brain.nii.gz')
    
    print("\n coregistering: "+brain_r1+"\n to:"+mni_t1)
    if not os.path.exists(brain_r1) or not os.path.exists(brain_r2s) or not os.path.exists(brain_qsm): 
        print("error: whole brain doesn't exist")
    elif not os.path.exists(mni_t1) or not os.path.exists(mni_t2) or not os.path.exists(mni_pd): 
        print("error: MNI template doesn't exist")
    else:
        mni_results = nr.registration.embedded_antsreg_multi(
                            source_images=[brain_r1,brain_r2s,brain_qsm],
                            target_images=[mni_t1,mni_pd,mni_t2],
                            run_rigid=True, run_affine=True,run_syn=True,
                            rigid_iterations=10000,
                            affine_iterations=2000,
                            coarse_iterations=240, 
                            medium_iterations=80, fine_iterations=40,
                            cost_function='MutualInformation', 
                            interpolation='NearestNeighbor',
                            regularization='Medium',
                            save_data=True, file_name='sub-'+subject+"_ses-1_acq-wb_mod-r1map_orient-std_brain_map-mni09b_reg-medium",
                            output_dir='/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/mni2009b-hr-mappings/',
                            overwrite=True)
