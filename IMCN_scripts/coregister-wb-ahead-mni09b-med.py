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
ahead_r1 = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/ahead_mni09b_med_r1map_n104.nii.gz'
ahead_r2s = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/ahead_mni09b_med_r2smap_n104.nii.gz'
ahead_qsm = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/ahead_mni09b_med_qsm_n104.nii.gz'

# Overall approach: 1) use nifti header + rigid registration to align whole brain and slab data, 
# 2) register whole brain to MNI (rigidly), 3) transfer slab data to MNI space. Coregistration 
# usually works better from whole brain to slab, so we do that and take the inverse mapping. 
# Additional labels for the whole brain data can also be mapped back onto the slab that way.

for subject in range(0,110): 
    subject = str(subject).zfill(3)
    wb_r1 = os.path.join('/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/r1maps/',
                            'sub-'+subject+'_ses-1_acq-wb_mod-r1map_orient-std_brain.nii.gz')
    wb_r2s = os.path.join('/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/r2smap/',
                            'sub-'+subject+'_ses-1_acq-wb_mod-r2starmap_orient-std_brain.nii.gz')
    wb_qsm = os.path.join('/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/qsm/',
                            'sub-'+subject+'_ses-1_acq-wb_mod-qsm_orient-std_brain.nii.gz')
    
    print("\n coregistering: "+wb_r1+"\n to:"+ahead_r1)
    print("\n coregistering: "+wb_r2s+"\n to:"+ahead_r2s)
    print("\n coregistering: "+wb_qsm+"\n to:"+ahead_qsm)
    if os.path.exists(wb_r1) is False: print("error: whole brain doesn't exist")
    elif os.path.exists(wb_r2s) is False: print("error: whole brain doesn't exist")
    elif os.path.exists(wb_qsm) is False: print("error: whole brain doesn't exist")
    elif os.path.exists(ahead_r1) is False: print("error: ahead template doesn't exist")
    elif os.path.exists(ahead_r2s) is False: print("error: ahead template doesn't exist")
    elif os.path.exists(ahead_qsm) is False: print("error: ahead template doesn't exist")
    else:
        mni_results = nr.registration.embedded_antsreg_multi(
                            source_images=[wb_r1,wb_r2s,wb_qsm],
                            target_images=[ahead_r1,ahead_r2s,ahead_qsm],
                            run_rigid=True, run_affine=True,run_syn=True,
                            rigid_iterations=10000,
                            affine_iterations=2000,
                            coarse_iterations=180, 
                            medium_iterations=60, fine_iterations=30,
                            cost_function='MutualInformation', 
                            interpolation='NearestNeighbor',
                            regularization='High',
                            save_data=True, file_name='sub-'+subject+"_ses-1_acq-wb_orient-std_brain_map-aheadMni09bMed",
                            output_dir='/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/mappings2aheadMni09bMed/',
                            overwrite=False)
