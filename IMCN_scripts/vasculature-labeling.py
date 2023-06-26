import nighres
import numpy
import nibabel
from glob import glob
import os
import subprocess
import scipy.ndimage

main_dir='/home/pilou/Projects/Ahead-Database/Vascular-Modeling/'
in_dir='/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/'
vessel_dir=main_dir+'process_qmri2/'
label_dir=main_dir+'vessel-labels/'

subjects = range(0,1)

for num in subjects:
    subject = 'sub-'+str(num).zfill(3)

    print("process subject "+subject)
    
    r1_file = in_dir+subject+'/ses-1/anat/wb/qmri/'+subject+'_ses-1_acq-wb_mod-r1hz_orient-std_brain.nii.gz'
    r2_file = in_dir+subject+'/ses-1/anat/wb/qmri/'+subject+'_ses-1_acq-wb_mod-r2hz_orient-std_brain.nii.gz'
    qsm_file = in_dir+subject+'/ses-1/anat/wb/qmri/'+subject+'_ses-1_acq-wb_mod-qsm_orient-std_brain.nii.gz'

    if (os.path.isfile(r1_file) and os.path.isfile(r2_file) and os.path.isfile(qsm_file)):

        print("Mask qR1")
        r1 = nighres.io.load_volume(r1_file)
        mask = scipy.ndimage.binary_erosion((r1.get_fdata()>0), iterations=5)
        mask_img = nibabel.Nifti1Image(mask, r1.affine, r1.header)
        
        print("1. Vasculature from R2*")
        
        vascR2 = nighres.filtering.multiscale_vessel_filter(r2_file, structure_intensity='bright', filterType = 'RRF',
                                        propagationtype = 'diffusion', threshold = 0.5,
                                        factor = 0.5,max_diff = 0.001,max_itr = 100,
                                        scale_step = 1.0,
                                        scales = 3,
                                        prior_image=mask_img,
                                        invert_prior=False,
                                        save_data=True,
                                        overwrite=False,
                                        output_dir=vessel_dir)   
        
        print("2. Vasculature from R1")
        
        vascR1 = nighres.filtering.multiscale_vessel_filter(r1_file, structure_intensity='bright', filterType = 'RRF',
                                        propagationtype = 'diffusion', threshold = 0.5,
                                        factor = 0.5,max_diff = 0.001,max_itr = 100,
                                        scale_step = 1.0,
                                        scales = 3,
                                        prior_image=mask_img,
                                        invert_prior=False,
                                        save_data=True,
                                        overwrite=False,
                                        output_dir=vessel_dir)   
        
        print("3. QSM")
            
        vascQSM = nighres.filtering.multiscale_vessel_filter(qsm_file, structure_intensity='bright', filterType = 'RRF',
                                        propagationtype = 'diffusion', threshold = 0.5,
                                        factor = 0.5,max_diff = 0.001,max_itr = 100,
                                        scale_step = 1.0,
                                        scales = 3,
                                        prior_image=mask_img,
                                        invert_prior=False,
                                        save_data=True,
                                        overwrite=False,
                                        output_dir=vessel_dir)   

        print("4. From voxels to segments")
        
        # seems a reasonable set of parameters, good groupings but still small segments
        clsR1 = nighres.microscopy.directional_line_clustering(labels=[vascR1['segmentation']],
                                                     directions=[vascR1['direction']],
                                                    scales=[vascR1['diameter']],
                                                    thickness=0.0,angle=60.0,distance=6.0,probability=0.5,anisotropy=0.5,
                                                    voxels=True,relabel=False,mip=0,
                                                    save_data=True,overwrite=False,
                                                    output_dir=label_dir,
                                                    file_name=subject+'_art-labels-lines')

        clsR2 = nighres.microscopy.directional_line_clustering(labels=[vascR2['segmentation']],
                                                     directions=[vascR2['direction']],
                                                    scales=[vascR2['diameter']],
                                                    thickness=0.0,angle=60.0,distance=6.0,probability=0.5,anisotropy=0.5,
                                                    voxels=True,relabel=False,mip=0,
                                                    save_data=True,overwrite=False,
                                                    output_dir=label_dir,
                                                    file_name=subject+'_vasc-labels-lines')

        clsQSM = nighres.microscopy.directional_line_clustering(labels=[vascQSM['segmentation']],
                                                     directions=[vascQSM['direction']],
                                                    scales=[vascQSM['diameter']],
                                                    thickness=0.0,angle=60.0,distance=6.0,probability=0.5,anisotropy=0.5,
                                                    voxels=True,relabel=False,mip=0,
                                                    save_data=True,overwrite=False,
                                                    output_dir=label_dir,
                                                    file_name=subject+'_ven-labels-lines')

        # global groupings?
        # 45 deg / 4 vox too exclusive, 60 deg / 6 vox too inclusive...
        # still some left-out regions due to the first stage of processing
        clsall = nighres.microscopy.directional_line_clustering(labels=[clsR1['grouping'],clsR2['grouping'],clsQSM['grouping']],
                                                    directions=[vascR1['direction'],vascR2['direction'],vascQSM['direction']],
                                                    scales=[vascR1['diameter'],vascR2['diameter'],vascQSM['diameter']],
                                                    thickness=0.0,angle=90.0,distance=6.0,probability=0.5,anisotropy=0.5,
                                                    voxels=False,relabel=False,across=True,mip=10,
                                                    save_data=True,overwrite=True,
                                                    output_dir=label_dir,
                                                    file_name=subject+'_vasc-labels-merged')

        # global voxel groupings?
        clsall = nighres.microscopy.directional_line_clustering(labels=[vascR1['segmentation'],vascR2['segmentation'],vascQSM['segmentation']],
                                                    directions=[vascR1['direction'],vascR2['direction'],vascQSM['direction']],
                                                    scales=[vascR1['diameter'],vascR2['diameter'],vascQSM['diameter']],
                                                    thickness=0.0,angle=90.0,distance=2.0,probability=0.5,anisotropy=0.5,
                                                    voxels=True,relabel=False,across=False,mip=10,
                                                    save_data=True,overwrite=True,
                                                    output_dir=label_dir,
                                                    file_name=subject+'_vasc-voxels-merged')
        