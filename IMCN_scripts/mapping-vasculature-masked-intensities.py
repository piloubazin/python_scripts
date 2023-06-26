import nighres
import numpy
import nibabel
from glob import glob
import os
import subprocess
import scipy.ndimage


in_dir='/home/pilou/Projects/Ahead-Database/Vascular-Modeling/process_qmri2/'
img_dir='/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/'
sub_dir='/home/pilou/Projects/Ahead-Database/Automated-Parcellation/qmri2/'
mask_dir='/home/pilou/Projects/Ahead-Database/Vascular-Modeling/BG_Artefact/'

out_dir='/home/pilou/Projects/Ahead-Database/Vascular-Modeling/process_mni09b/'
tmp_dir='/home/pilou/Projects/Ahead-Database/Vascular-Modeling/tmp/'

# change to get a different sub-set, can also use a dictionary of IDs
subjects = range(0,110)
#subjects = [0,1,2]

# statistics
for num in subjects:
    subject = 'sub-'+str(num).zfill(3)

    print("process subject "+subject)
    
    # here we are taking the partial volume of vessels detected on the R2*
    # and the corresponding subcortical structures
    # (other combinations are possible)
    vasc_seg_file = in_dir+subject+'_ses-1_acq-wb_mod-r2hz_orient-std_brain_mvf-seg.nii.gz'
    vasc_dia_file = in_dir+subject+'_ses-1_acq-wb_mod-r2hz_orient-std_brain_mvf-dia.nii.gz'
    vasc_pv_file = in_dir+subject+'_ses-1_acq-wb_mod-r2hz_orient-std_brain_mvf-pv.nii.gz'
    
    massp_seg_file = sub_dir+subject+'_atlas10_31struct_qmri2_massp-label.nii.gz'
    
    r2s_img_file = img_dir+subject+'/ses-1/anat/wb/qmri/'+subject+'_ses-1_acq-wb_mod-r2hz_orient-std_brain.nii.gz'

    artefact_file = glob(mask_dir+subject+'*.nii.gz')
    

    if (os.path.isfile(vasc_seg_file) and os.path.isfile(r2s_img_file) and os.path.isfile(massp_seg_file)):

        # mask the artifact values
        mask = nighres.io.load_volume(artefact_file[0]).get_fdata()
        
        dia = nighres.io.load_volume(vasc_dia_file)
        dia = nibabel.Nifti1Image((1.0-mask)*dia.get_fdata(), dia.affine, dia.header)
        
        seg = nighres.io.load_volume(vasc_seg_file)
        seg = nibabel.Nifti1Image((1.0-mask)*seg.get_fdata(), seg.affine, seg.header)
        
        pv = nighres.io.load_volume(vasc_pv_file)
        pv = nibabel.Nifti1Image((1.0-mask)*pv.get_fdata(), pv.affine, pv.header)
        
        
        
        # intensity stats
        nighres.statistics.segmentation_statistics(segmentation=massp_seg_file, intensity=None, template=None,
                                statistics=["Volume"], 
                                output_csv='vascular-seg-r2s-noartefact-volume-statistics.csv', file_name=subject+'_stats',
                                atlas=None, skip_first=True, ignore_zero=False)
        
        nighres.statistics.segmentation_statistics(segmentation=massp_seg_file, intensity=seg, template=None,
                                statistics=["Sum_intensity"], 
                                output_csv='vascular-seg-r2s-noartefact-vessel-length-statistics.csv', file_name=subject+'_stats',
                                atlas=None, skip_first=True, ignore_zero=True)
                
        nighres.statistics.segmentation_statistics(segmentation=massp_seg_file, intensity=dia, template=None,
                                statistics=["Mean_intensity", "Std_intensity", "Median_intensity", "IQR_intensity"], 
                                output_csv='vascular-seg-r2s-noartefact-diameter-statistics.csv', file_name=subject+'_stats',
                                atlas=None, skip_first=True, ignore_zero=True)
                
        nighres.statistics.segmentation_statistics(segmentation=massp_seg_file, intensity=pv, template=None,
                                statistics=["Sum_intensity", "Mean_intensity" "Median_intensity"], 
                                output_csv='vascular-seg-r2s-noartefact-pv-statistics.csv', file_name=subject+'_stats',
                                atlas=None, skip_first=True, ignore_zero=False)
                
        nighres.statistics.segmentation_statistics(segmentation=massp_seg_file, intensity=r2s_img_file, template=None,
                                statistics=["Mean_intensity", "Std_intensity", "Median_intensity", "IQR_intensity"], 
                                output_csv='vascular-seg-r2s-noartefact-r2s-statistics.csv', file_name=subject+'_stats',
                                atlas=None, skip_first=True, ignore_zero=False)
                
        