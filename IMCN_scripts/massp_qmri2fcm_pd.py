import nighres
import numpy
import nibabel
from glob import glob
import os
import subprocess
import math
import scipy.ndimage

in_dir = '/home/pilou/Projects/PD-Ahead/qMRI-Recomputed-good-skullstrip/'
out_dir = '/home/pilou/Projects/PD-Ahead/Massp-PD/massp-qmri2fcm/'

subjects = ['sub-111','sub-112','sub-113','sub-114','sub-115','sub-116','sub-117','sub-118','sub-119','sub-120','sub-121']

for subject in subjects:

    proc_dir=out_dir+'processing/'+subject+'/'
    sub_dir=out_dir+subject+'/subcortex/'
    ctx_dir=out_dir+subject+'/cortex/'
    vsl_dir=out_dir+subject+'/vessels/'
    tis_dir=out_dir+subject+'/tissues/'
    
    # input names
    qr1_file = in_dir+subject+'/ses-1/anat/wb/qmri/'+subject+'_ses-1_acq-wb2_mod-r1hz_orient-std_brain.nii.gz'
    qr2s_file = in_dir+subject+'/ses-1/anat/wb/qmri/'+subject+'_ses-1_acq-wb2_mod-r2hz_orient-std_brain.nii.gz'
    qsm_file = in_dir+subject+'/ses-1/anat/wb/qmri/'+subject+'_ses-1_acq-wb2_mod-qsm_orient-std_brain.nii.gz'
    qpd_file = in_dir+subject+'/ses-1/anat/wb/qmri/'+subject+'_ses-1_acq-wb2_mod-qpd_orient-std_brain.nii.gz'
    
    # output suffixes, limited to what is strictly necessary
    sub_label = '_massp-labels.nii.gz'
    sub_proba = '_massp-probas.nii.gz'
    
    mgdm_label = '_mgdm-labels.nii.gz'
    mgdm_proba = '_mgdm-probas.nii.gz'
    
    ctx_label = '_cruise-cortex.nii.gz'
    ctx_gwb = '_cruise-wmsurf.nii.gz'
    ctx_cgb = '_cruise-pialsurf.nii.gz'
    
    vsl_pv = '_mvf-pv.nii.gz'
    vsl_dia = '_mvf-dia.nii.gz'
    
    tis_class = '_rfcm-class.nii.gz'
    tis_map = '_rfcm-map.nii.gz'
    
    # main processing steps
    
    print("\n0. Tissue intensity maps\n")
    
    if (not os.path.exists(tis_dir)):
        os.makedirs(tis_dir)
    
    tis_map_file = tis_dir+subject+'_ses-1_r1tim'+tis_map
    if (not os.path.exists(tis_map_file)):
    
        rfcm = nighres.segmentation.fuzzy_cmeans(qr1_file, clusters=2, max_iterations=150, max_difference=0.001, 
                        smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                        save_data=True, overwrite=False, output_dir=proc_dir)
    
        os.rename(rfcm['intensity'],tis_map_file)
    qr1_file = tis_map_file
    
    tis_map_file = tis_dir+subject+'_ses-1_r2stim'+tis_map
    if (not os.path.exists(tis_map_file)):
    
        rfcm = nighres.segmentation.fuzzy_cmeans(qr2s_file, clusters=2, max_iterations=150, max_difference=0.001, 
                        smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                        save_data=True, overwrite=False, output_dir=proc_dir)
    
        os.rename(rfcm['intensity'],tis_map_file)
    qr2s_file = tis_map_file
    
    tis_map_file = tis_dir+subject+'_ses-1_qsmtim'+tis_map
    if (not os.path.exists(tis_map_file)):
    
        rfcm = nighres.segmentation.fuzzy_cmeans(qsm_file, clusters=3, max_iterations=150, max_difference=0.001, 
                        smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                        save_data=True, overwrite=False, output_dir=proc_dir)
    
        os.rename(rfcm['intensity'],tis_map_file)
    qsm_file = tis_map_file
    
    tis_map_file = tis_dir+subject+'_ses-1_qpdtim'+tis_map
    if (not os.path.exists(tis_map_file)):
    
        rfcm = nighres.segmentation.fuzzy_cmeans(qpd_file, clusters=3, max_iterations=150, max_difference=0.001, 
                        smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                        save_data=True, overwrite=False, output_dir=proc_dir)
    
        os.rename(rfcm['intensity'],tis_map_file)
    qpd_file = tis_map_file
    
    print("\n1. MASSP subcortical parcellation\n")
    
    if (not os.path.exists(sub_dir)):
        os.makedirs(sub_dir)
    
    sub_label_file = sub_dir+subject+'_ses-1'+sub_label
    sub_proba_file = sub_dir+subject+'_ses-1'+sub_proba
    
    if (not os.path.exists(sub_label_file) or not os.path.exists(sub_proba_file)):
        template = nighres.data.download_AHEAD_template()
        
#        ants = nighres.registration.embedded_antsreg_multi(
#                                source_images=[qr1_file,qr2s_file,qsm_file],
#                                target_images=[template['qr1'],template['qr2s'],template['qsm']],
#                                run_rigid=True, run_affine=True, run_syn=True,
#                                rigid_iterations=10000,
#                                affine_iterations=2000,
#                                coarse_iterations=180, 
#                                medium_iterations=60, fine_iterations=30,
#                                cost_function='MutualInformation', 
#                                interpolation='Linear',
#                                regularization='Medium',
#                                ignore_affine=True, 
#                                save_data=True, file_name=subject+'_ses-1_reg2ahead',
#                                output_dir=proc_dir)
        
        # two-step ANTs registration for maximum precision
        ants1 = nighres.registration.embedded_antsreg_multi(
                                    source_images=[qr1_file,qr2s_file],
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
                                    save_data=True, file_name=subject+'_ses-1_reg2ahead-step1',
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
                                    save_data=True, file_name=subject+'_ses-1_reg2ahead-step2',
                                    output_dir=proc_dir)
    
        mapping = nighres.registration.apply_coordinate_mappings(ants1['mapping'],mapping1=ants2['mapping'],
                                    save_data=True, overwrite=False, file_name=subject+'_ses-1_reg2ahead-mapping',
                                    output_dir=proc_dir)
    
        inverse = nighres.registration.apply_coordinate_mappings(ants2['inverse'],mapping1=ants1['inverse'],
                                    save_data=True, overwrite=False, file_name=subject+'_ses-1_reg2ahead-inverse',
                                    output_dir=proc_dir)
    

        # not sure which option is best...
    #    massp_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/massp_atlases/'
    #    hist = nighres.segmentation.conditional_shape_map_intensities(31, 3, 3,
    #                          contrast_images=[template['qr1'],template['qr2s'],template['qsm']], 
    #                          target_images=[ants['transformed_sources'][0],ants['transformed_sources'][1],ants['transformed_sources'][2]],
    #                          shape_atlas_probas=massp_dir+'massp_17structures_spatial_proba.nii.gz', 
    #                          shape_atlas_labels=massp_dir+'massp_17structures_spatial_label.nii.gz',
    #                          skeleton_atlas_probas=massp_dir+'massp_17structures_skeleton_proba.nii.gz', 
    #                          skeleton_atlas_labels=massp_dir+'massp_17structures_skeleton_proba.nii.gz',
    #                          intensity_atlas_hist=massp_dir+'massp_17structures_qmri2_r1r2sqsm_histograms.nii.gz',
    #                          save_data=True, output_dir=proc_dir, file_name=subject+'_massp_hist')
        
    #    massp = nighres.parcellation.massp(target_images=[qr1_file,qr2s_file,qsm_file],
    #                                    map_to_target=ants['inverse'], 
    #                                    #intensity_atlas_hist='/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/massp_atlases/massp_17structures_qmri2_r1r2sqsm_histograms.nii.gz',
    #                                    intensity_atlas_hist=hist['cond_hist'],
    #                                    max_iterations=120, max_difference=0.1,
    #                                    save_data=True, file_name=subject+'_ses-1_reg2ahead',
    #                                    output_dir=proc_dir, overwrite=False)
    
        atlas_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/qmri2_ahead/'
        massp = nighres.parcellation.massp(target_images=[qr1_file,qr2s_file,qsm_file,qpd_file],
                                    map_to_target=inverse['result'], 
                                    shape_atlas_probas=atlas_dir+'atlas10_31struct_qmri2fcm_massp-sproba.nii.gz', 
                                    shape_atlas_labels=atlas_dir+'atlas10_31struct_qmri2fcm_massp-slabel.nii.gz',
                                    intensity_atlas_hist=atlas_dir+'atlas10_31struct_qmri2fcm_massp-chist.nii.gz',
                                    skeleton_atlas_probas=atlas_dir+'atlas10_31struct_qmri2fcm_massp-kproba.nii.gz', 
                                    skeleton_atlas_labels=atlas_dir+'atlas10_31struct_qmri2fcm_massp-klabel.nii.gz',
                                    max_iterations=120, max_difference=0.1,
                                    #atlas_file='/home/pilou/Code/github/nighres/nighres/parcellation/massp_17structures_manual.json',
                                    intensity_prior=0.25,
                                    save_data=True, file_name=subject+'_ses-1_tim',
                                    output_dir=proc_dir, overwrite=False)
    
        os.rename(massp['max_label'],sub_label_file)
        os.rename(massp['max_proba'],sub_proba_file)
    
    
    print("\n2. MGDM-CRUISE pipeline for cortex\n")
    
    if (not os.path.exists(ctx_dir)):
        os.makedirs(ctx_dir)
    
    mgdm_label_file = ctx_dir+subject+'_ses-1'+mgdm_label
    mgdm_proba_file = ctx_dir+subject+'_ses-1'+mgdm_proba
    
    if (not os.path.exists(mgdm_label_file) or not os.path.exists(mgdm_proba_file)):
        qt1_file = proc_dir+subject+'_ses-1_acq-mp2rage_mod-t1msnob1_orient-std_brain.nii.gz'
        if (not os.path.isfile(qt1_file)):
            r1 = nighres.io.load_volume(qr1_file)
            data = r1.get_fdata()
            data = 1000.0*numpy.divide(1.0,data,where=data>0)
            img = nibabel.Nifti1Image(data,r1.affine,r1.header)
            nighres.io.save_volume(qt1_file, img)
        
        mgdm = nighres.brain.mgdm_segmentation(contrast_image1=qr1_file,
                                                contrast_type1="R1FCM",
                                                atlas_file="brain-atlas-quant-3.0.9.txt",
                                                adjust_intensity_priors=False,
                                                normalize_qmaps=False,
                                                compute_posterior=False, posterior_scale=5.0,
                                                diffuse_probabilities=True,
                                                save_data=True, file_name=subject+'_ses-1',
                                                output_dir=proc_dir)
    
        os.rename(mgdm['segmentation'],mgdm_label_file)
        os.rename(mgdm['distance'],mgdm_proba_file)
        
        mgdm['segmentation'] = mgdm_label_file
        mgdm['distance'] = mgdm_proba_file
    else:
        mgdm_mems_file = proc_dir+subject+'_ses-1_mgdm-mems.nii.gz'
        mgdm_lbls_file = proc_dir+subject+'_ses-1_mgdm-lbls.nii.gz'
        
        mgdm = {'segmentation': mgdm_label_file,
                'distance': mgdm_proba_file,
                'memberships': mgdm_mems_file,
                'labels': mgdm_lbls_file}
    
    regions = ['left_cerebrum','right_cerebrum','cerebellum']
    rgns = ['lcr','rcr','cb']
    for idx,region in enumerate(regions):
        ctx_label_file = ctx_dir+subject+'_ses-1_'+rgns[idx]+ctx_label
        ctx_gwb_file = ctx_dir+subject+'_ses-1_'+rgns[idx]+ctx_gwb
        ctx_cgb_file = ctx_dir+subject+'_ses-1_'+rgns[idx]+ctx_cgb
    
        if (not os.path.exists(ctx_label_file) or not os.path.exists(ctx_gwb_file) or not os.path.exists(ctx_cgb_file)):
    
            cortex = nighres.brain.extract_brain_region(segmentation=mgdm['segmentation'],
                                                        levelset_boundary=mgdm['distance'],
                                                        maximum_membership=mgdm['memberships'],
                                                        maximum_label=mgdm['labels'],
                                                        extracted_region=region,
                                                        save_data=True,
                                                        output_dir=proc_dir, 
                                                        file_name=subject+'_ses-1_'+rgns[idx])
        
            cruise = nighres.cortex.cruise_cortex_extraction(
                                    init_image=cortex['inside_mask'],
                                    wm_image=cortex['inside_proba'],
                                    gm_image=cortex['region_proba'],
                                    csf_image=cortex['background_proba'],
                                    normalize_probabilities=True,
                                    save_data=True, output_dir=proc_dir,
                                    file_name=subject+'_ses-1_'+rgns[idx])
    
            os.rename(cruise['cortex'],ctx_label_file)
            os.rename(cruise['gwb'],ctx_gwb_file)
            os.rename(cruise['cgb'],ctx_cgb_file)
             
            
    
    print("\n3. Extraction of the vasculature\n")
    
    if (not os.path.exists(vsl_dir)):
        os.makedirs(vsl_dir)
    
    vsl_pv_file = vsl_dir+subject+'_ses-1_r2vsl'+vsl_pv
    vsl_dia_file = vsl_dir+subject+'_ses-1_r2vsl'+vsl_dia
    
    art_pv_file = vsl_dir+subject+'_ses-1_r1art'+vsl_pv
    art_dia_file = vsl_dir+subject+'_ses-1_r1art'+vsl_dia
    
    ven_pv_file = vsl_dir+subject+'_ses-1_qsmven'+vsl_pv
    ven_dia_file = vsl_dir+subject+'_ses-1_qsmven'+vsl_dia
    
    if (not os.path.exists(vsl_pv_file) or not os.path.exists(vsl_dia_file) \
        or not os.path.exists(art_pv_file) or not os.path.exists(art_dia_file) \
        or not os.path.exists(ven_pv_file) or not os.path.exists(ven_dia_file)) :
        
        # we erode the mask to remove the false (and true) positives at the border (dura, etc)
        r1 = nighres.io.load_volume(qr1_file)
        mask = scipy.ndimage.binary_erosion((r1.get_fdata()>0), iterations=5)
        mask_img = nibabel.Nifti1Image(mask, r1.affine, r1.header)
                
        vsl = nighres.filtering.multiscale_vessel_filter(qr2s_file, structure_intensity='bright', filterType = 'RRF',
                                        propagationtype = 'diffusion', threshold = 0.5,
                                        factor = 0.5,max_diff = 0.001,max_itr = 100,
                                        scale_step = 1.0,
                                        scales = 3,
                                        prior_image=mask_img,
                                        invert_prior=False,
                                        save_data=True,
                                        overwrite=False,
                                        output_dir=proc_dir)   
                
        art = nighres.filtering.multiscale_vessel_filter(qr1_file, structure_intensity='bright', filterType = 'RRF',
                                        propagationtype = 'diffusion', threshold = 0.5,
                                        factor = 0.5,max_diff = 0.001,max_itr = 100,
                                        scale_step = 1.0,
                                        scales = 3,
                                        prior_image=mask_img,
                                        invert_prior=False,
                                        save_data=True,
                                        overwrite=False,
                                        output_dir=proc_dir)   
                
        ven = nighres.filtering.multiscale_vessel_filter(qsm_file, structure_intensity='bright', filterType = 'RRF',
                                        propagationtype = 'diffusion', threshold = 0.5,
                                        factor = 0.5,max_diff = 0.001,max_itr = 100,
                                        scale_step = 1.0,
                                        scales = 3,
                                        prior_image=mask_img,
                                        invert_prior=False,
                                        save_data=True,
                                        overwrite=False,
                                        output_dir=proc_dir)   
        
        os.rename(vsl['pv'],vsl_pv_file)
        os.rename(vsl['diameter'],vsl_dia_file)
    
        os.rename(art['pv'],art_pv_file)
        os.rename(art['diameter'],art_dia_file)
    
        os.rename(ven['pv'],ven_pv_file)
        os.rename(ven['diameter'],ven_dia_file)
    
