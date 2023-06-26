import nibabel as nb
import nighres as nr
import glob
import numpy as np
import os
from nighres.io import load_volume, save_volume


ahead_r1 = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/ahead_qmri2_mni09b_med_r1map_n105.nii.gz'
ahead_r2s = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/ahead_qmri2_mni09b_med_r2map_n105.nii.gz'
ahead_qsm = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/ahead_qmri2_mni09b_med_qsmap_n105.nii.gz'
ahead_qpd = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/ahead_qmri2_mni09b_med_pdmap_n105.nii.gz'

# Overall approach: 1) use nifti header + rigid registration to align whole brain and slab data, 
# 2) register whole brain to MNI (rigidly), 3) transfer slab data to MNI space. Coregistration 
# usually works better from whole brain to slab, so we do that and take the inverse mapping. 
# Additional labels for the whole brain data can also be mapped back onto the slab that way.

for subject in range(0,110): 
    subject = str(subject).zfill(3)
    wb_r1 = '/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/sub-'+subject\
    			+'/ses-1/anat/wb/qmri/sub-'+subject+'_ses-1_acq-wb_mod-r1hz_orient-std_brain.nii.gz'
    wb_r2s = '/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/sub-'+subject\
    			+'/ses-1/anat/wb/qmri/sub-'+subject+'_ses-1_acq-wb_mod-r2hz_orient-std_brain.nii.gz'
    wb_qsm = '/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/sub-'+subject\
    			+'/ses-1/anat/wb/qmri/sub-'+subject+'_ses-1_acq-wb_mod-qsm_orient-std_brain.nii.gz'
    wb_qpd = '/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/sub-'+subject\
    			+'/ses-1/anat/wb/qmri/sub-'+subject+'_ses-1_acq-wb_mod-qpd_orient-std_brain.nii.gz'
    
    print("\n coregistering: "+wb_r1+"\n to:"+ahead_r1)
    print("\n coregistering: "+wb_r2s+"\n to:"+ahead_r2s)
    print("\n coregistering: "+wb_qsm+"\n to:"+ahead_qsm)
    print("\n coregistering: "+wb_qpd+"\n to:"+ahead_qpd)
    if os.path.exists(wb_r1) is False: print("error: whole brain r1 doesn't exist")
    elif os.path.exists(wb_r2s) is False: print("error: whole brain r2* doesn't exist")
    elif os.path.exists(wb_qsm) is False: print("error: whole brain qsm doesn't exist")
    elif os.path.exists(wb_qpd) is False: print("error: whole brain qpd doesn't exist")
    elif os.path.exists(ahead_r1) is False: print("error: ahead template r1 doesn't exist")
    elif os.path.exists(ahead_r2s) is False: print("error: ahead template r2* doesn't exist")
    elif os.path.exists(ahead_qsm) is False: print("error: ahead template qsm doesn't exist")
    elif os.path.exists(ahead_qpd) is False: print("error: ahead template qpd doesn't exist")
    else:
        mni_results = nr.registration.embedded_antsreg_multi(
                            source_images=[wb_r1,wb_r2s,wb_qsm,wb_qpd],
                            target_images=[ahead_r1,ahead_r2s,ahead_qsm,ahead_qpd],
                            run_rigid=True, run_affine=True,run_syn=True,
                            rigid_iterations=10000,
                            affine_iterations=2000,
                            coarse_iterations=180, 
                            medium_iterations=60, fine_iterations=30,
                            cost_function='MutualInformation', 
                            interpolation='NearestNeighbor',
                            regularization='High',
                            save_data=True, file_name='sub-'+subject+"_ses-1_acq-wb_orient-std_brain_map-ahead-qmri2",
                            output_dir='/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/ahead-qmri2-mappings/',
                            overwrite=False)
