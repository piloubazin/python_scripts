import nighres
import numpy
import nibabel
from glob import glob
import os
import subprocess
import math
import scipy.ndimage

# hopefully we just have to change the subject id every time
in_dir = '/home/pilou/Projects/IronSleep/Mpm7T-Pipeline/data/'
out_dir = '/home/pilou/Projects/IronSleep/Mpm7T-Pipeline/derivatives/nighres/'
proc_dir = '/home/pilou/Projects/IronSleep/Mpm7T-Pipeline/processing/'

subjects = ['sub-li2am220707115619n0001','sub-li2am220707115619n0002',
            'sub-li2am220707115619n0003','sub-li2am220707115619n0004']
session = 'hMRI_maps/'
thr_dir='threshold/'
tis_dir='tissues/'
sub_dir='subcortex/'
ctx_dir='cortex/'
vsl_dir='vessels/'

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

thr_img = '_thr.nii.gz'


contrast = {'MT':'_MPM_MTsat_masked','PD':'_MPM_PD_masked',
            'R1':'_MPM_R1_masked','R2s':'_MPM_R2s_OLS_masked'}
            
for subject in subjects:

    # clean up unusual values
    if (not os.path.exists(out_dir+subject+'/'+thr_dir)):
        os.makedirs(out_dir+subject+'/'+thr_dir)

    # note: keep zero as lower bound for masking
    r1th_file = out_dir+subject+'/'+thr_dir+subject+contrast['R1']+thr_img
    if not os.path.isfile(r1th_file):
        print("Threshold R1 map")
        r1 = nighres.io.load_volume(in_dir+subject+'/'+session+subject+contrast['R1']+'.nii.gz')
        r1th = nibabel.Nifti1Image(numpy.minimum(2.0,numpy.maximum(0.0,r1.get_fdata())), r1.affine, r1.header)
        r1th = nibabel.as_closest_canonical(r1th)
        nighres.io.save_volume(r1th_file, r1th)
    
    r2th_file = out_dir+subject+'/'+thr_dir+subject+contrast['R2s']+thr_img
    if not os.path.isfile(r2th_file):
        print("Threshold R2s map")
        r2 = nighres.io.load_volume(in_dir+subject+'/'+session+subject+contrast['R2s']+'.nii.gz')
        r2th = nibabel.Nifti1Image(numpy.minimum(200.0,numpy.maximum(0.0,r2.get_fdata())), r2.affine, r2.header)
        r2th = nibabel.as_closest_canonical(r2th)
        nighres.io.save_volume(r2th_file, r2th)
    
    pdth_file = out_dir+subject+'/'+thr_dir+subject+contrast['PD']+thr_img
    if not os.path.isfile(pdth_file):
        print("Threshold PD map")
        pd = nighres.io.load_volume(in_dir+subject+'/'+session+subject+contrast['PD']+'.nii.gz')
        pdth = nibabel.Nifti1Image(numpy.minimum(200.0,numpy.maximum(0.0,pd.get_fdata())), pd.affine, pd.header)
        pdth = nibabel.as_closest_canonical(pdth)
        nighres.io.save_volume(pdth_file, pdth)  
    
    mtth_file = out_dir+subject+'/'+thr_dir+subject+contrast['MT']+thr_img
    if not os.path.isfile(mtth_file):
        print("Threshold MT map")
        mt = nighres.io.load_volume(in_dir+subject+'/'+session+subject+contrast['MT']+'.nii.gz')
        mtth = nibabel.Nifti1Image(numpy.minimum(2.0,numpy.maximum(0.0,mt.get_fdata())), mt.affine, mt.header)
        mtth = nibabel.as_closest_canonical(mtth)
        nighres.io.save_volume(mtth_file, mtth)
    

    # main processing steps
    
    print("\n0. Tissue intensity maps\n")
    
    if (not os.path.exists(out_dir+subject+'/'+tis_dir)):
        os.makedirs(out_dir+subject+'/'+tis_dir)
    
    tis_map_file = out_dir+subject+'/'+tis_dir+subject+'_r1tim'+tis_map
    if (not os.path.exists(tis_map_file)):
    
        rfcm = nighres.segmentation.fuzzy_cmeans(r1th_file, clusters=2, max_iterations=150, max_difference=0.001, 
                        smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                        save_data=True, overwrite=False, output_dir=proc_dir)
    
        os.rename(rfcm['intensity'],tis_map_file)
    qr1_file = tis_map_file
    
    tis_map_file = out_dir+subject+'/'+tis_dir+subject+'_r2stim'+tis_map
    if (not os.path.exists(tis_map_file)):
    
        rfcm = nighres.segmentation.fuzzy_cmeans(r2th_file, clusters=2, max_iterations=150, max_difference=0.001, 
                        smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                        save_data=True, overwrite=False, output_dir=proc_dir)
    
        os.rename(rfcm['intensity'],tis_map_file)
    qr2s_file = tis_map_file
    
    tis_map_file = out_dir+subject+'/'+tis_dir+subject+'_mttim'+tis_map
    if (not os.path.exists(tis_map_file)):
    
        rfcm = nighres.segmentation.fuzzy_cmeans(mtth_file, clusters=2, max_iterations=150, max_difference=0.001, 
                        smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                        save_data=True, overwrite=False, output_dir=proc_dir)
    
        os.rename(rfcm['intensity'],tis_map_file)
    qmt_file = tis_map_file
    
    tis_map_file = out_dir+subject+'/'+tis_dir+subject+'_pdtim'+tis_map
    if (not os.path.exists(tis_map_file)):
    
        rfcm = nighres.segmentation.fuzzy_cmeans(pdth_file, clusters=3, max_iterations=150, max_difference=0.001, 
                        smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                        save_data=True, overwrite=False, output_dir=proc_dir)
    
        os.rename(rfcm['intensity'],tis_map_file)
    qpd_file = tis_map_file
    
    print("\n1. MASSP subcortical parcellation\n")
    
    if (not os.path.exists(out_dir+subject+'/'+sub_dir)):
        os.makedirs(out_dir+subject+'/'+sub_dir)
    
    sub_label_file = out_dir+subject+'/'+sub_dir+subject+sub_label
    sub_proba_file = out_dir+subject+'/'+sub_dir+subject+sub_proba
    
    if (not os.path.exists(sub_label_file) or not os.path.exists(sub_proba_file)):
        template = ['/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/ahead_qmri2_mni09b_med_r1map_n105.nii.gz',
                    '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/ahead_qmri2_mni09b_med_r2map_n105.nii.gz',
                    '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/ahead_qmri2_mni09b_med_pdmap_n105.nii.gz']
        
        # two-step ANTs registration for maximum precision
        ants1 = nighres.registration.embedded_antspy_multi(
                                    source_images=[qr1_file,qr2s_file,qpd_file],
                                    target_images=[template[0],template[1],template[2]],
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
                                    save_data=True, file_name=subject+'_reg2ahead-step1',
                                    output_dir=proc_dir)
            
        ants2 = nighres.registration.embedded_antspy_multi(
                                    source_images=[ants1['transformed_sources'][0],ants1['transformed_sources'][1],ants1['transformed_sources'][2]],
                                    target_images=[template[0],template[1],template[2]],
                                    run_rigid=False, run_affine=False, run_syn=True,
                                    coarse_iterations=180, 
                                    medium_iterations=60, fine_iterations=30,
                                    cost_function='MutualInformation', 
                                    interpolation='Linear',
                                    regularization='Medium',
                                    mask_zero=False,
                                    smooth_mask=0.33,
                                    ignore_affine=True, 
                                    save_data=True, file_name=subject+'_reg2ahead-step2',
                                    output_dir=proc_dir)
    
        mapping = nighres.registration.apply_coordinate_mappings(ants1['mapping'],mapping1=ants2['mapping'],
                                    save_data=True, overwrite=False, file_name=subject+'_reg2ahead-mapping',
                                    output_dir=proc_dir)
    
        inverse = nighres.registration.apply_coordinate_mappings(ants2['inverse'],mapping1=ants1['inverse'],
                                    save_data=True, overwrite=False, file_name=subject+'_reg2ahead-inverse',
                                    output_dir=proc_dir)
        
        # transform other contrasts for more elaborate mapping? or just R1, R2* for now?
    

        # mapping of contrast to MPMs, or use direclty fcm-normalized ones?
        #massp_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/massp_atlases/'
        #hist = nighres.segmentation.conditional_shape_map_intensities(31, 3, 3,
        #                      contrast_images=[template['qr1'],template['qr2s'],template['qpd']], 
        #                      target_images=[ants['transformed_sources'][0],ants['transformed_sources'][1],ants['transformed_sources'][2]],
        #                      shape_atlas_probas=massp_dir+'massp_17structures_spatial_proba.nii.gz', 
        #                      shape_atlas_labels=massp_dir+'massp_17structures_spatial_label.nii.gz',
        #                      skeleton_atlas_probas=massp_dir+'massp_17structures_skeleton_proba.nii.gz', 
        #                      skeleton_atlas_labels=massp_dir+'massp_17structures_skeleton_proba.nii.gz',
        #                      intensity_atlas_hist=massp_dir+'massp_17structures_qmri2_r1r2sqsm_histograms.nii.gz',
        #                      save_data=True, output_dir=proc_dir, file_name=subject+'_massp_mpm_hist')
        
        #massp = nighres.parcellation.massp(target_images=[qr1_file,qr2s_file,qpd_file],
        #                                map_to_target=ants['inverse'], 
        #                                intensity_atlas_hist=hist['cond_hist'],
        #                                max_iterations=120, max_difference=0.1,
        #                                save_data=True, file_name=subject+'_reg2ahead',
        #                                output_dir=proc_dir, overwrite=False)
    
        atlas_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/qmri2_ahead/'
        massp = nighres.parcellation.massp(target_images=[qr1_file,qr2s_file,qpd_file],
                                    map_to_target=inverse['result'], 
                                    shape_atlas_probas=atlas_dir+'atlas10_31struct_qmri2fcm_massp-sproba.nii.gz', 
                                    shape_atlas_labels=atlas_dir+'atlas10_31struct_qmri2fcm_massp-slabel.nii.gz',
                                    intensity_atlas_hist=atlas_dir+'atlas10_31struct_qmri2fcm_massp-chist-r1-r2s-qpd.nii.gz',
                                    skeleton_atlas_probas=atlas_dir+'atlas10_31struct_qmri2fcm_massp-kproba.nii.gz', 
                                    skeleton_atlas_labels=atlas_dir+'atlas10_31struct_qmri2fcm_massp-klabel.nii.gz',
                                    max_iterations=120, max_difference=0.1,
                                    intensity_prior=0.25,
                                    save_data=True, file_name=subject+'_tim',
                                    output_dir=proc_dir, overwrite=False)
    
        os.rename(massp['max_label'],sub_label_file)
        os.rename(massp['max_proba'],sub_proba_file)
    
    
    print("\n2. MGDM-CRUISE pipeline for cortex\n")
    
    if (not os.path.exists(out_dir+subject+'/'+ctx_dir)):
        os.makedirs(out_dir+subject+'/'+ctx_dir)
    
    mgdm_label_file = out_dir+subject+'/'+ctx_dir+subject+mgdm_label
    mgdm_proba_file = out_dir+subject+'/'+ctx_dir+subject+mgdm_proba
    
    if (not os.path.exists(mgdm_label_file) or not os.path.exists(mgdm_proba_file)):
        qt1_file = out_dir+subject+'/'+thr_dir+subject+contrast['R1']+thr_img
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
                                                save_data=True, file_name=subject,
                                                output_dir=proc_dir)
    
        os.rename(mgdm['segmentation'],mgdm_label_file)
        os.rename(mgdm['distance'],mgdm_proba_file)
        
        mgdm['segmentation'] = mgdm_label_file
        mgdm['distance'] = mgdm_proba_file
    else:
        mgdm_mems_file = proc_dir+subject+'_mgdm-mems.nii.gz'
        mgdm_lbls_file = proc_dir+subject+'_mgdm-lbls.nii.gz'
        
        mgdm = {'segmentation': mgdm_label_file,
                'distance': mgdm_proba_file,
                'memberships': mgdm_mems_file,
                'labels': mgdm_lbls_file}
    
    regions = ['left_cerebrum','right_cerebrum','cerebellum']
    rgns = ['lcr','rcr','cb']
    for idx,region in enumerate(regions):
        ctx_label_file = out_dir+subject+'/'+ctx_dir+subject+rgns[idx]+ctx_label
        ctx_gwb_file = out_dir+subject+'/'+ctx_dir+subject+rgns[idx]+ctx_gwb
        ctx_cgb_file = out_dir+subject+'/'+ctx_dir+subject+rgns[idx]+ctx_cgb
    
        if (not os.path.exists(ctx_label_file) or not os.path.exists(ctx_gwb_file) or not os.path.exists(ctx_cgb_file)):
    
            cortex = nighres.brain.extract_brain_region(segmentation=mgdm['segmentation'],
                                                        levelset_boundary=mgdm['distance'],
                                                        maximum_membership=mgdm['memberships'],
                                                        maximum_label=mgdm['labels'],
                                                        extracted_region=region,
                                                        save_data=True,
                                                        output_dir=proc_dir, 
                                                        file_name=subject+'_'+rgns[idx])
        
            cruise = nighres.cortex.cruise_cortex_extraction(
                                    init_image=cortex['inside_mask'],
                                    wm_image=cortex['inside_proba'],
                                    gm_image=cortex['region_proba'],
                                    csf_image=cortex['background_proba'],
                                    normalize_probabilities=True,
                                    save_data=True, output_dir=proc_dir,
                                    file_name=subject+'_'+rgns[idx])
    
            os.rename(cruise['cortex'],ctx_label_file)
            os.rename(cruise['gwb'],ctx_gwb_file)
            os.rename(cruise['cgb'],ctx_cgb_file)
             
            
    
    print("\n3. Extraction of the vasculature\n")
    
    if (not os.path.exists(out_dir+subject+'/'+vsl_dir)):
        os.makedirs(out_dir+subject+'/'+vsl_dir)
    
    vsl_pv_file = out_dir+subject+'/'+vsl_dir+subject+'_r2vsl'+vsl_pv
    vsl_dia_file = out_dir+subject+'/'+vsl_dir+subject+'_r2vsl'+vsl_dia
    
    art_pv_file = out_dir+subject+'/'+vsl_dir+subject+'_r1art'+vsl_pv
    art_dia_file = out_dir+subject+'/'+vsl_dir+subject+'_r1art'+vsl_dia
    
    if (not os.path.exists(vsl_pv_file) or not os.path.exists(vsl_dia_file) \
        or not os.path.exists(art_pv_file) or not os.path.exists(art_dia_file)) :
        
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
                
        os.rename(vsl['pv'],vsl_pv_file)
        os.rename(vsl['diameter'],vsl_dia_file)
    
        os.rename(art['pv'],art_pv_file)
        os.rename(art['diameter'],art_dia_file)
