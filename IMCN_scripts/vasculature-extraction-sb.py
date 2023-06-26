import nighres
import numpy
import nibabel
from glob import glob
import os
import subprocess
import scipy.ndimage

main_dir='/home/pilou/Projects/Ahead-Database/Vascular-Modeling/'
in_dir='/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/'
out_dir=main_dir+'process_qmri2sb/'
final_dir=main_dir+'vessel-maps/'

subjects = range(91,92)

for num in subjects:
    subject = 'sub-'+str(num).zfill(3)

    print("process subject "+subject)
    
    r1_file = in_dir+subject+'/ses-1/anat/sb/qmri/'+subject+'_ses-1_acq-sb2_mod-r1hz_orient-std_brain.nii.gz'
    r2_file = in_dir+subject+'/ses-1/anat/sb/qmri/'+subject+'_ses-1_acq-sb2_mod-r2hz_orient-std_brain.nii.gz'
    qsm_file = in_dir+subject+'/ses-1/anat/sb/qmri/'+subject+'_ses-1_acq-sb2_mod-qsm_orient-std_brain.nii.gz'

    if (os.path.isfile(r1_file) and os.path.isfile(r2_file) and os.path.isfile(qsm_file)):

        print("Mask qR1")
        r1 = nighres.io.load_volume(r1_file)
        mask = scipy.ndimage.binary_erosion((r1.get_fdata()>0), iterations=5)
        mask_img = nibabel.Nifti1Image(mask, r1.affine, r1.header)
        
        print("1. Vasculature from R2*")
        
        vasc = nighres.filtering.multiscale_vessel_filter(r2_file, structure_intensity='bright', filterType = 'RRF',
                                        propagationtype = 'diffusion', threshold = 0.5,
                                        factor = 0.5,max_diff = 0.001,max_itr = 100,
                                        scale_step = 1.0,
                                        scales = 3,
                                        prior_image=mask_img,
                                        invert_prior=False,
                                        save_data=True,
                                        overwrite=False,
                                        output_dir=out_dir)   
        
        print("2. Vasculature from R1")
        
        vasc = nighres.filtering.multiscale_vessel_filter(r1_file, structure_intensity='bright', filterType = 'RRF',
                                        propagationtype = 'diffusion', threshold = 0.5,
                                        factor = 0.5,max_diff = 0.001,max_itr = 100,
                                        scale_step = 1.0,
                                        scales = 3,
                                        prior_image=mask_img,
                                        invert_prior=False,
                                        save_data=True,
                                        overwrite=False,
                                        output_dir=out_dir)   
        
        print("3. QSM")
            
        vasc = nighres.filtering.multiscale_vessel_filter(qsm_file, structure_intensity='bright', filterType = 'RRF',
                                        propagationtype = 'diffusion', threshold = 0.5,
                                        factor = 0.5,max_diff = 0.001,max_itr = 100,
                                        scale_step = 1.0,
                                        scales = 3,
                                        prior_image=mask_img,
                                        invert_prior=False,
                                        save_data=True,
                                        overwrite=False,
                                        output_dir=out_dir)   
