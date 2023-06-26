import nighres
import os
import nibabel
from glob import glob
import numpy
import scipy.ndimage

## use previously computed example as testbed
atlas_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/qmri2_ahead/'
data_dir='/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/'
proc_dir ='/home/pilou/Projects/Ahead-Database/Automated-Parcellation/qmri2fcm_proc/'
out_dir ='/home/pilou/Projects/Ahead-Database/Automated-Parcellation/qmri2fcm_all/'


#for subject in atlas_subjects: 
#for subject_number in range(0,110): 
for subject_number in range(70,71): 
    subject = 'sub-'+str(subject_number).zfill(3)
    
    session = 'ses-2'

    r1_subject = glob(data_dir+subject+'/'+session+'/anat/wb/qmri/'+subject+'*wb2*r1hz_orient-std_brain.nii.gz')
    r2s_subject = glob(data_dir+subject+'/'+session+'/anat/wb/qmri/'+subject+'*wb2*r2hz_orient-std_brain.nii.gz')
    qsm_subject = glob(data_dir+subject+'/'+session+'/anat/wb/qmri/'+subject+'*wb2*qsm_orient-std_brain.nii.gz')
    qpd_subject = glob(data_dir+subject+'/'+session+'/anat/wb/qmri/'+subject+'*wb2*qpd_orient-std_brain.nii.gz')
    
    ## check for data availability ##
    if len(r1_subject)>0 and len(r2s_subject)>0 and len(qsm_subject)>0 and len(qpd_subject)>0:
        
        r1_subject = r1_subject[0]
        r2s_subject = r2s_subject[0]
        qsm_subject = qsm_subject[0]
        qpd_subject = qpd_subject[0]
        
        # need to erode qsm mask
        qsm_img = nighres.io.load_volume(qsm_subject)
        brainmask = (qsm_img.get_fdata()!=0)
        qsmmask = scipy.ndimage.binary_erosion(brainmask, iterations=5)
        qsm_img = nibabel.Nifti1Image(qsmmask*qsm_img.get_fdata(), qsm_img.affine, qsm_img.header)
        nighres.io.save_volume(proc_dir+subject+'_'+session+'_mod-qsm_masked.nii.gz', qsm_img)
        qsm_img = proc_dir+subject+'_'+session+'_mod-qsm_masked.nii.gz'
            
        # generate fcm-based contrasts
        rfcm = nighres.segmentation.fuzzy_cmeans(r1_subject, clusters=2, max_iterations=150, max_difference=0.001, 
                                                 smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                 save_data=True, overwrite=False, output_dir=proc_dir)
        r1_subject = rfcm['intensity']
        
        rfcm = nighres.segmentation.fuzzy_cmeans(r2s_subject, clusters=2, max_iterations=150, max_difference=0.001, 
                                                 smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                 save_data=True, overwrite=False, output_dir=proc_dir)
        r2s_subject = rfcm['intensity']
        
        rfcm = nighres.segmentation.fuzzy_cmeans(qsm_img, clusters=3, max_iterations=150, max_difference=0.001, 
                                                 smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                 save_data=True, overwrite=False, output_dir=proc_dir)
        qsm_subject = rfcm['intensity']
        
        rfcm = nighres.segmentation.fuzzy_cmeans(qpd_subject, clusters=3, max_iterations=150, max_difference=0.001, 
                                                 smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                 save_data=True, overwrite=False, output_dir=proc_dir)
        qpd_subject = rfcm['intensity']
        
        
        # two-step ANTs registration for maximum precision
        template = nighres.data.download_AHEAD_template()

        ants1 = nighres.registration.embedded_antsreg_multi(
                                    source_images=[r1_subject,r2s_subject],
                                    target_images=[template['qr1'],template['qr2s']],
                                    run_rigid=True, run_affine=True, run_syn=True,
                                    rigid_iterations=10000,
                                    affine_iterations=2000,
                                    coarse_iterations=180, 
                                    medium_iterations=60, fine_iterations=30,
                                    cost_function='MutualInformation', 
                                    interpolation='Linear',
                                    regularization='High',
                                    smooth_mask=0.33,
                                    ignore_affine=True, 
                                    save_data=True, file_name=subject+'_'+session+'_reg2ahead-step1',
                                    output_dir=proc_dir)
            
        ants2 = nighres.registration.embedded_antsreg_multi(
                            source_images=[ants1['transformed_sources'][0],ants1['transformed_sources'][1]],
                            target_images=[template['qr1'],template['qr2s']],
                            run_rigid=False, run_affine=False, run_syn=True,
                            coarse_iterations=180, 
                            medium_iterations=60, fine_iterations=30,
                            cost_function='MutualInformation', 
                            interpolation='Linear',
                            regularization='Medium',
                            mask_zero=False,
                            smooth_mask=0.33,
                            ignore_affine=True, 
                            save_data=True, file_name=subject+'_'+session+'_reg2ahead-step2',
                            output_dir=proc_dir)
    
        mapping = nighres.registration.apply_coordinate_mappings(ants1['mapping'],mapping1=ants2['mapping'],
                                    save_data=True, overwrite=False, file_name=subject+'_'+session+'_reg2ahead-mapping',
                                    output_dir=proc_dir)
    
        inverse = nighres.registration.apply_coordinate_mappings(ants2['inverse'],mapping1=ants1['inverse'],
                                    save_data=True, overwrite=False, file_name=subject+'_'+session+'_reg2ahead-inverse',
                                    output_dir=proc_dir)
    
        massp = nighres.parcellation.massp(target_images=[r1_subject,r2s_subject,qsm_subject,qpd_subject],
                                    map_to_target=inverse['result'], 
                                    shape_atlas_probas=atlas_dir+'atlas10_31struct_qmri2fcm_massp-sproba.nii.gz', 
                                    shape_atlas_labels=atlas_dir+'atlas10_31struct_qmri2fcm_massp-slabel.nii.gz',
                                    intensity_atlas_hist=atlas_dir+'atlas10_31struct_qmri2fcm_massp-chist.nii.gz',
                                    skeleton_atlas_probas=atlas_dir+'atlas10_31struct_qmri2fcm_massp-kproba.nii.gz', 
                                    skeleton_atlas_labels=atlas_dir+'atlas10_31struct_qmri2fcm_massp-klabel.nii.gz',
                                    max_iterations=120, max_difference=0.1,
                                    #atlas_file='/home/pilou/Code/github/nighres/nighres/parcellation/massp_17structures_manual.json',
                                    intensity_prior=0.25,
                                    save_data=True, file_name=subject+'_'+session+'_tim',
                                    output_dir=proc_dir, overwrite=False)

 