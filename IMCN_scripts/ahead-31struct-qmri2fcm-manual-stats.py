import nighres
import os
import nibabel
import numpy
from glob import glob

## use previously computed example as testbed
data_dir='/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/'
lbl_dir ='/home/pilou/Projects/Ahead-Database/Automated-Parcellation/manual_delineations/'
thk_dir ='/home/pilou/Projects/Ahead-Database/Automated-Parcellation/manual_delineations/thickness/'

big5 = ['stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r']        

for target_number in range(0,110): 
    target = 'sub-'+str(target_number).zfill(3)

    r1_target = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*_acq-wb2_mod-r1hz_orient-std_brain.nii.gz')
    r2s_target = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*acq-wb2_mod-r2hz_orient-std_brain.nii.gz')
    qsm_target = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*acq-wb2_mod-qsm_orient-std_brain.nii.gz')
    #qpd_target = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*acq-wb2_mod-qpd_orient-std_brain.nii.gz')

    ## check for data availability ##
    if len(r1_target)>0 and len(r2s_target)>0 and len(qsm_target)>0:
        
        r1_target = r1_target[0]
        r2s_target = r2s_target[0]
        qsm_target = qsm_target[0]
        #qpd_target = qpd_target[0]
        
        for lbl in big5:
            mask = glob(lbl_dir+target+'*_mask-'+lbl+'_rat*.nii.gz')
            if len(mask)>0: 
                mask = mask[0]
                mask_img = nighres.io.load_volume(mask)
                mask_img = nibabel.Nifti1Image(numpy.flip(mask_img.get_fdata(),axis=0),mask_img.affine,mask_img.header)
                
                nighres.statistics.segmentation_statistics(segmentation=mask_img, intensity=r1_target, template=None,
                                statistics=["Median_intensity"], 
                                output_csv='ahead-'+lbl+'_qmri2fcm-manual_r1hz-med.csv', file_name=target+'_stats',
                                atlas=None, skip_first=True, ignore_zero=True)
                
                nighres.statistics.segmentation_statistics(segmentation=mask_img, intensity=r2s_target, template=None,
                                statistics=["Median_intensity"], 
                                output_csv='ahead-'+lbl+'_qmri2fcm-manual_r2hz-med.csv', file_name=target+'_stats',
                                atlas=None, skip_first=True, ignore_zero=True)
                
                nighres.statistics.segmentation_statistics(segmentation=mask_img, intensity=qsm_target, template=None,
                                statistics=["Median_intensity"], 
                                output_csv='ahead-'+lbl+'_qmri2fcm-manual_qsm-med.csv', file_name=target+'_stats',
                                atlas=None, skip_first=True, ignore_zero=True)
                
                thickness = nighres.shape.levelset_thickness(mask_img, shape_image_type='probability_map',
                                save_data=False,overwrite=False,output_dir=thk_dir)['thickness']
                
                nighres.statistics.segmentation_statistics(segmentation=mask_img, intensity=thickness, template=None,
                                statistics=["Median_intensity"], 
                                output_csv='ahead-'+lbl+'_qmri2fcm-manual_thickness-med.csv', file_name=target+'_stats',
                                atlas=None, skip_first=True, ignore_zero=True)

                nighres.statistics.segmentation_statistics(segmentation=mask_img, intensity=None, template=None,
                                statistics=["Volume"], 
                                output_csv='ahead-'+lbl+'_qmri2fcm-manual_vol.csv', file_name=target+'_stats',
                                atlas=None, skip_first=True, ignore_zero=True)

