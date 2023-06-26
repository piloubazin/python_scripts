import nighres
import nibabel
import os
import glob
import numpy

data_dir='/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/'
out_dir = '/home/pilou/Projects/Ahead-Database/IronMyelin-mapping/'

for target_number in range(0,110): 
    target = 'sub-'+str(target_number).zfill(3)

    r1_target = glob.glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*_acq-wb2_mod-r1hz_orient-std_brain.nii.gz')
    r2s_target = glob.glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*_acq-wb2_mod-r2hz_orient-std_brain.nii.gz')
    qsm_target = glob.glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*_acq-wb2_mod-qsm_orient-std_brain.nii.gz')
    qpd_target = glob.glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*_acq-wb2_mod-qpd_orient-std_brain.nii.gz')
    
    ## check for data availability ##
    if len(r1_target)>0:
        
        r1_target = r1_target[0]
        r2s_target = r2s_target[0]
        qsm_target = qsm_target[0]
        qpd_target = qpd_target[0]
        
        header = nighres.io.load_volume(r1_target).header
        affine = nighres.io.load_volume(r1_target).affine
        r1_data = nighres.io.load_volume(r1_target).get_fdata()
        r2s_data = nighres.io.load_volume(r2s_target).get_fdata()
        qpd_data = nighres.io.load_volume(qpd_target).get_fdata()
        qsm_data = nighres.io.load_volume(qsm_target).get_fdata()
        
        # plug in the model
        # (old version)
        #iron_data = -2.0935*qpd_data + 142.2635*qsm_data + 0.2148*r2s_data + 3.1705*r1_data
        #myelin_data = -3.2042*qpd_data + -109.2097*qsm_data + -0.0487*r2s_data + 28.8318*r1_data
        # new version
        #iron_data = 5.3789 + 125.5717*qsm_data + 0.2633*r2s_data
        #myelin_data = -7.1401 + 35.7341*r1_data - 0.1919*r2s_data
        # latest version (wb2)
        #iron_data = -3.4565 + 126.1282*qsm_data + 0.2328*r2s_data
        #myelin_data = -8.3468 + 36.3114*r1_data - 0.1850*r2s_data
        # improved version (with medians, etc)
        #iron_data = -3.690616 + 119.094908*qsm_data + 0.240898*r2s_data
        #myelin_data = -7.728833 + 31.920727*r1_data - 0.121142*r2s_data
        # improved version with tighter registration (qmri2fcm)
        iron_data = -3.828341 + 98.279472*qsm_data + 0.243101*r2s_data
        myelin_data = -7.989650 + 31.874833*r1_data + -0.111816*r2s_data

        iron_img = nibabel.Nifti1Image(iron_data,affine,header)
        nighres.io.save_volume(out_dir+target+'_qmri2map-iron.nii.gz', iron_img)
        
        myelin_img = nibabel.Nifti1Image(myelin_data,affine,header)
        nighres.io.save_volume(out_dir+target+'_qmri2map-myelin.nii.gz', myelin_img)
        