import nighres
import os
import nibabel
from glob import glob
import numpy
import scipy.ndimage

## use previously computed example as testbed
data_dir='/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/'
in_dir ='/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/'
out_dir = '/home/pilou/Projects/Ahead-Database/Slab-Analysis/qmri2fcm_slab10/'
atlas_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/qmri2fcm_ahead/'
proc_dir ='/home/pilou/Projects/Ahead-Database/Automated-Parcellation/qmri2fcm_proc/'


subjects = ['sub-000','sub-001','sub-002','sub-004','sub-008','sub-012','sub-024','sub-031','sub-033','sub-040']
targets = ['sub-000','sub-001','sub-002','sub-004','sub-008','sub-012','sub-024','sub-031','sub-033','sub-040']
#targets = ['sub-000']

wb_structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r',
              'ic_hem-l','ic_hem-r','vta_hem-l','vta_hem-r','fx_hem-lr','ppn_hem-l','ppn_hem-r','cl_hem-l','cl_hem-r','lh_hem-l','lh_hem-r',
              'ico_hem-l','ico_hem-r','pag_hem-l','pag_hem-r','sco_hem-l','sco_hem-r','scg_hem-l','scg_hem-r']
wb_structures = []

sb_structures = ['ac_hem-lr','pc_hem-lr','chn_hem-l','chn_hem-r',
              'lgn_hem-l','lgn_hem-r','mgn_hem-l','mgn_hem-r','ns_hem-l','ns_hem-r']

session = 'ses-1'

for target in targets:

    n_contrasts = 4
    n_subjects = len(subjects)-1
    n_structures = len(wb_structures)+len(sb_structures)

    r1_target = glob(data_dir+target+'/ses-1/anat/sb/qmri/'+target+'*sb2*r1hz_orient-std_brain.nii.gz')[0]
    r2s_target = glob(data_dir+target+'/ses-1/anat/sb/qmri/'+target+'*sb2*r2hz_orient-std_brain.nii.gz')[0]
    qsm_target = glob(data_dir+target+'/ses-1/anat/sb/qmri/'+target+'*sb2*qsm_orient-std_brain.nii.gz')[0]
    qpd_target = glob(data_dir+target+'/ses-1/anat/sb/qmri/'+target+'*sb2*qpd_orient-std_brain.nii.gz')[0]
    
    r1_wb = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*r1hz_orient-std_brain.nii.gz')[0]
    r2s_wb = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*r2hz_orient-std_brain.nii.gz')[0]
    qsm_wb = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*qsm_orient-std_brain.nii.gz')[0]
    qpd_wb = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*qpd_orient-std_brain.nii.gz')[0]
    
    # need to erode qsm mask??
    #qsm_img = nighres.io.load_volume(qsm_target)
    #brainmask = (qsm_img.get_fdata()!=0)
    #qsmmask = scipy.ndimage.binary_erosion(brainmask, iterations=5)
    #qsm_img = nibabel.Nifti1Image(qsmmask*qsm_img.get_fdata(), qsm_img.affine, qsm_img.header)
    #nighres.io.save_volume(atlas_dir+target+'_mod-qsm_masked.nii.gz', qsm_img)
    #qsm_img = atlas_dir+target+'_mod-qsm_masked.nii.gz'
        
    # generate fcm-based contrasts: maybe not a good idea here?
    #rfcm = nighres.segmentation.fuzzy_cmeans(r1_target, clusters=2, max_iterations=150, max_difference=0.001, 
    #                                         smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
    #                                         save_data=True, overwrite=False, output_dir=atlas_dir)
    
    rfcm = nighres.segmentation.fuzzy_cmeans(r1_target, clusters=3, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir, 
                                             file_name=target+'_mod-r1hz_clusters-3.nii.gz')
    r1_target = rfcm['intensity']
    
    #rfcm = nighres.segmentation.fuzzy_cmeans(r2s_target, clusters=2, max_iterations=150, max_difference=0.001, 
    #                                         smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
    #                                         save_data=True, overwrite=False, output_dir=atlas_dir)
    
    rfcm = nighres.segmentation.fuzzy_cmeans(r2s_target, clusters=3, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir, 
                                             file_name=target+'_mod-r2hz_clusters-3.nii.gz')
    r2s_target = rfcm['intensity']
    
    rfcm = nighres.segmentation.fuzzy_cmeans(qsm_target, clusters=3, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir)
    qsm_target = rfcm['intensity']
    
    rfcm = nighres.segmentation.fuzzy_cmeans(qpd_target, clusters=3, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir)
    qpd_target = rfcm['intensity']

    # same for the whole brain version, for consistency
    # need to erode qsm mask??
    #qsm_img = nighres.io.load_volume(qsm_wb)
    #brainmask = (qsm_img.get_fdata()!=0)
    #qsmmask = scipy.ndimage.binary_erosion(brainmask, iterations=5)
    #qsm_img = nibabel.Nifti1Image(qsmmask*qsm_img.get_fdata(), qsm_img.affine, qsm_img.header)
    #nighres.io.save_volume(atlas_dir+target+'_mod-qsm_masked.nii.gz', qsm_img)
    #qsm_img = atlas_dir+target+'_mod-qsm_masked.nii.gz'
    
    #rfcm = nighres.segmentation.fuzzy_cmeans(r1_wb, clusters=2, max_iterations=150, max_difference=0.001, 
    #                                         smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
    #                                         save_data=True, overwrite=False, output_dir=atlas_dir)
    #r1_wb = rfcm['intensity']
    
    #rfcm = nighres.segmentation.fuzzy_cmeans(r2s_wb, clusters=2, max_iterations=150, max_difference=0.001, 
    #                                         smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
    #                                         save_data=True, overwrite=False, output_dir=atlas_dir)
    #r2s_wb = rfcm['intensity']
    
    #rfcm = nighres.segmentation.fuzzy_cmeans(qsm_img, clusters=3, max_iterations=150, max_difference=0.001, 
    #                                         smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
    #                                         save_data=True, overwrite=False, output_dir=atlas_dir)
    #qsm_wb = rfcm['intensity']
    
    #rfcm = nighres.segmentation.fuzzy_cmeans(qpd_wb, clusters=3, max_iterations=150, max_difference=0.001, 
    #                                         smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
    #                                         save_data=True, overwrite=False, output_dir=atlas_dir)
    #qpd_wb = rfcm['intensity']

    # explicitly recompute a better coreg
    template = nighres.data.download_AHEAD_template()
    
    # first computing the coreg to whole brain
    ants0 = nighres.registration.embedded_antspy_multi(
                                    source_images=[r1_target,r2s_target],
                                    target_images=[r1_wb,r2s_wb],
                                    run_rigid=True, run_affine=False, run_syn=False,
                                    rigid_iterations=10000,
                                    affine_iterations=2000,
                                    coarse_iterations=180, 
                                    medium_iterations=60, fine_iterations=30,
                                    cost_function='MutualInformation', 
                                    interpolation='Linear',
                                    regularization='High',
                                    smooth_mask=0.33,
                                    ignore_affine=False, 
                                    save_data=True, file_name=target+'_'+session+'_reg2slab-step0',
                                    output_dir=proc_dir)

    ants1 = nighres.registration.embedded_antspy_multi(
                                    source_images=[r1_wb,r2s_wb],
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
                                    save_data=True, file_name=target+'_'+session+'_reg2ahead-step1',
                                    output_dir=proc_dir)
            
    ants2 = nighres.registration.embedded_antspy_multi(
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
                        save_data=True, file_name=target+'_'+session+'_reg2ahead-step2',
                        output_dir=proc_dir)

    mapping = nighres.registration.apply_coordinate_mappings(ants0['mapping'],mapping1=ants1['mapping'],mapping2=ants2['mapping'],
                                save_data=True, overwrite=False, file_name=target+'_'+session+'_slab2ahead-mapping',
                                output_dir=proc_dir)

    inverse = nighres.registration.apply_coordinate_mappings(ants2['inverse'],mapping1=ants1['inverse'],mapping2=ants0['inverse'],
                                save_data=True, overwrite=False, file_name=target+'_'+session+'_slab2ahead-inverse',
                                output_dir=proc_dir)

    map_target = mapping['result']
    invmap_target = inverse['result']
    
    manual_target = glob(in_dir+'label_maps/'+target+'_labeling-slab10_v2.nii.gz')[0]
    
    target_contrasts = [r1_target,r2s_target,qsm_target,qpd_target]
    
    atlas_contrasts = []
    atlas_levelsets = []
    atlas_skeletons = []
    atlas_mappings = []
    for i in range(n_subjects):
        atlas_contrasts.append([])
        atlas_levelsets.append([])
        atlas_skeletons.append([])
        atlas_mappings.append([])
    
    i=0
    for sub in subjects:
        if sub is not target:
            # find the files
            r1_atlas = glob(data_dir+sub+'/ses-1/anat/sb/qmri/'+sub+'*sb2*r1hz_orient-std_brain.nii.gz')[0]
            r2s_atlas = glob(data_dir+sub+'/ses-1/anat/sb/qmri/'+sub+'*sb2*r2hz_orient-std_brain.nii.gz')[0]
            qsm_atlas = glob(data_dir+sub+'/ses-1/anat/sb/qmri/'+sub+'*sb2*qsm_orient-std_brain.nii.gz')[0]
            qpd_atlas = glob(data_dir+sub+'/ses-1/anat/sb/qmri/'+sub+'*sb2*qpd_orient-std_brain.nii.gz')[0]
            
            r1_awb = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2*r1hz_orient-std_brain.nii.gz')[0]
            r2s_awb = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2*r2hz_orient-std_brain.nii.gz')[0]
            qsm_awb = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2*qsm_orient-std_brain.nii.gz')[0]
            qpd_awb = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2*qpd_orient-std_brain.nii.gz')[0]
            
            # need to erode qsm mask
            #qsm_img = nighres.io.load_volume(qsm_atlas)
            #brainmask = (qsm_img.get_fdata()!=0)
            #qsmmask = scipy.ndimage.binary_erosion(brainmask, iterations=5)
            #qsm_img = nibabel.Nifti1Image(qsmmask*qsm_img.get_fdata(), qsm_img.affine, qsm_img.header)
            #nighres.io.save_volume(atlas_dir+sub+'_mod-qsm_masked.nii.gz', qsm_img)
            #qsm_img = atlas_dir+sub+'_mod-qsm_masked.nii.gz'
            
            # generate fcm-based contrasts
            #rfcm = nighres.segmentation.fuzzy_cmeans(r1_atlas, clusters=2, max_iterations=150, max_difference=0.001, 
            #                                         smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
            #                                         save_data=True, overwrite=False, output_dir=atlas_dir)
            rfcm = nighres.segmentation.fuzzy_cmeans(r1_atlas, clusters=3, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir, 
                                                     file_name=sub+'_mod-r1hz_clusters-3.nii.gz')
            r1_atlas = rfcm['intensity']
            
            #rfcm = nighres.segmentation.fuzzy_cmeans(r2s_atlas, clusters=2, max_iterations=150, max_difference=0.001, 
            #                                         smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
            #                                         save_data=True, overwrite=False, output_dir=atlas_dir)
            rfcm = nighres.segmentation.fuzzy_cmeans(r2s_atlas, clusters=3, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir, 
                                                     file_name=sub+'_mod-r2hz_clusters-3.nii.gz')
            r2s_atlas = rfcm['intensity']
            
            rfcm = nighres.segmentation.fuzzy_cmeans(qsm_atlas, clusters=3, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir)
            qsm_atlas = rfcm['intensity']
            
            rfcm = nighres.segmentation.fuzzy_cmeans(qpd_atlas, clusters=3, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir)
            qpd_atlas = rfcm['intensity']
            
            ants0 = nighres.registration.embedded_antspy_multi(
                                    source_images=[r1_atlas,r2s_atlas],
                                    target_images=[r1_awb,r2s_awb],
                                    run_rigid=True, run_affine=False, run_syn=False,
                                    rigid_iterations=10000,
                                    affine_iterations=2000,
                                    coarse_iterations=180, 
                                    medium_iterations=60, fine_iterations=30,
                                    cost_function='MutualInformation', 
                                    interpolation='Linear',
                                    regularization='High',
                                    smooth_mask=0.33,
                                    ignore_affine=False, 
                                    save_data=True, file_name=sub+'_'+session+'_reg2slab-step0',
                                    output_dir=proc_dir)

            ants1 = nighres.registration.embedded_antspy_multi(
                                            source_images=[r1_awb,r2s_awb],
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
                                            save_data=True, file_name=sub+'_'+session+'_reg2ahead-step1',
                                            output_dir=proc_dir)
                    
            ants2 = nighres.registration.embedded_antspy_multi(
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
                                save_data=True, file_name=sub+'_'+session+'_reg2ahead-step2',
                                output_dir=proc_dir)
        
            mapping = nighres.registration.apply_coordinate_mappings(ants0['mapping'],mapping1=ants1['mapping'],mapping2=ants2['mapping'],
                                        save_data=True, overwrite=False, file_name=sub+'_'+session+'_slab2ahead-mapping',
                                        output_dir=proc_dir)
        
            inverse = nighres.registration.apply_coordinate_mappings(ants2['inverse'],mapping1=ants1['inverse'],mapping2=ants0['inverse'],
                                        save_data=True, overwrite=False, file_name=sub+'_'+session+'_slab2ahead-inverse',
                                        output_dir=proc_dir)
        
            mapping_wb = nighres.registration.apply_coordinate_mappings(ants1['mapping'],mapping1=ants2['mapping'],
                                        save_data=True, overwrite=False, file_name=sub+'_'+session+'_reg2ahead-mapping',
                                        output_dir=proc_dir)
        
            map_atlas = mapping['result']
            invmap_atlas = inverse['result']
            map_wba = mapping['result']
            
            print("atlas subject: "+sub)
        
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=r1_atlas,
                            mapping1=map_atlas, interpolation='linear',
                            #save_data=True, file_name=sub+'_mod-r1hz_def-qmri2fcm_slab-atlas.nii.gz',
                            save_data=True, file_name=sub+'_mod-r1hz_clusters-3_def-qmri2fcm_slab-atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
        
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=r2s_atlas,
                            mapping1=map_atlas, interpolation='linear',
                            #save_data=True, file_name=sub+'_mod-r2hz_def-qmri2fcm_slab-atlas.nii.gz',
                            save_data=True, file_name=sub+'_mod-r2hz_clusters-3_def-qmri2fcm_slab-atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
        
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=qsm_atlas,
                            mapping1=map_atlas, interpolation='linear',
                            save_data=True, file_name=sub+'_mod-qsm_def-qmri2fcm_slab-atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
            
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=qpd_atlas,
                            mapping1=map_atlas, interpolation='linear',
                            save_data=True, file_name=sub+'_mod-qpd_def-qmri2fcm_slab-atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
            
            atlas_mappings[i] = map_atlas
        
            for j,struct in enumerate(wb_structures):
                
                print('structure: '+struct+' (subject: '+sub+')')
                
                mask = glob(in_dir+'manual_delineations/best/'+sub+'*mask-'+struct+'*.nii.gz')[0]
                # flip the masks because of the change of canonical orientation
                # only to be done for old masks, check sources for new ones, if any
                mask_img = nighres.io.load_volume(mask)
                mask_img = nibabel.Nifti1Image(numpy.flip(mask_img.get_fdata(),axis=0),mask_img.affine,mask_img.header)

                defmask = nighres.registration.apply_coordinate_mappings(image=mask_img,
                                    mapping1=map_wba, interpolation='linear',
                                    save_data=True, file_name=sub+'_mask-'+struct+'_def-qmri2_atlas.nii.gz',
                                    output_dir=atlas_dir)['result']
                atlas_levelsets[i].append(nighres.surface.probability_to_levelset(
                                        probability_image=defmask,
                                        save_data=True, file_name=sub+'_mask-'+struct+'_def-qmri2_atlas.nii.gz',
                                        output_dir=atlas_dir)['result'])
                atlas_skeletons[i].append(nighres.shape.levelset_thickness(
                                        atlas_levelsets[i][j],
                                        save_data=True, file_name=sub+'_mask-'+struct+'_def-qmri2_atlas.nii.gz',
                                        output_dir=atlas_dir)['dist'])
                
            for j,struct in enumerate(sb_structures):
                
                print('structure: '+struct+' (subject: '+sub+')')
                
                mask = glob(in_dir+'manual_delineations/best/'+sub+'*mask-'+struct+'*.nii.gz')[0]
                # do NOT flip the masks that were done on the new slab data!!
                mask_img = nighres.io.load_volume(mask)
                #mask_img = nibabel.Nifti1Image(numpy.flip(mask_img.get_fdata(),axis=0),mask_img.affine,mask_img.header)
                
                # note: some slabs were oriented incorrectly in voxel space: align first
                mask_img = nibabel.as_closest_canonical(mask_img)

                defmask = nighres.registration.apply_coordinate_mappings(image=mask_img,
                                    mapping1=map_atlas, interpolation='linear',
                                    save_data=True, file_name=sub+'_mask-'+struct+'_def-qmri2_slab-atlas.nii.gz',
                                    output_dir=atlas_dir)['result']
                atlas_levelsets[i].append(nighres.surface.probability_to_levelset(
                                        probability_image=defmask,
                                        save_data=True, file_name=sub+'_mask-'+struct+'_def-qmri2_slab-atlas.nii.gz',
                                        output_dir=atlas_dir)['result'])
                atlas_skeletons[i].append(nighres.shape.levelset_thickness(
                                        atlas_levelsets[i][j],
                                        save_data=True, file_name=sub+'_mask-'+struct+'_def-qmri2_slab-atlas.nii.gz',
                                        output_dir=atlas_dir)['dist'])
                
            i = i+1;
    
    # this run is only to build the atlas
    nighres.parcellation.massp_atlasing(levelset_images=atlas_levelsets,
                            contrast_images=atlas_contrasts, 
                            skeleton_images=atlas_skeletons,
                            subjects=n_subjects, structures=n_structures,
                            contrasts=n_contrasts,
                            save_data=True, file_name='atlas91_slab_qmri2fcm_cl3_coregx2-'+target,
                            output_dir=atlas_dir, overwrite=False)
    
    nighres.parcellation.massp(target_images=target_contrasts,
                            map_to_target=invmap_target, structures=n_structures,
                            shape_atlas_probas=atlas_dir+'atlas91_slab_qmri2fcm_cl3_coregx2-'+target+'_massp-sproba.nii.gz', 
                            shape_atlas_labels=atlas_dir+'atlas91_slab_qmri2fcm_cl3_coregx2-'+target+'_massp-slabel.nii.gz',
                            intensity_atlas_hist=atlas_dir+'atlas91_slab_qmri2fcm_cl3_coregx2-'+target+'_massp-chist.nii.gz',
                            skeleton_atlas_probas=atlas_dir+'atlas91_slab_qmri2fcm_cl3_coregx2-'+target+'_massp-kproba.nii.gz', 
                            skeleton_atlas_labels=atlas_dir+'atlas91_slab_qmri2fcm_cl3_coregx2-'+target+'_massp-klabel.nii.gz',
                            max_iterations=120, max_difference=0.1,
                            intensity_prior=0.25,
                            save_data=True, file_name=target+'_atlas91_slab_qmri2fcm_cl3_coregx2',
                            output_dir=out_dir, overwrite=False)
    
    seg_result = out_dir+target+'_atlas91_slab_qmri2fcm_cl3_coregx2_massp-label.nii.gz'
    
    manual_img = nighres.io.load_volume(manual_target)
    #manual_img = nibabel.Nifti1Image(numpy.flip(manual_img.get_fdata(),axis=0),manual_img.affine,manual_img.header)
# skip
#    nighres.statistics.segmentation_statistics(segmentation=seg_result, intensity=None, template=manual_img,
#                                statistics=["Volumes","Volume_difference","Center_distance","Dice_overlap","Dilated_Dice_overlap","Average_surface_distance"], 
#                                output_csv='atlas91_slab10_qmri2fcm_coregx2_overlap_v2_cl3.csv', file_name=target,
#                                atlas=None, skip_first=True, ignore_zero=True)

    # test if results are different in atlas space
    manual_img = nighres.registration.apply_coordinate_mappings(
                            image=manual_img,
                            mapping1=map_target, interpolation='nearest',
                            save_data=True, file_name=target+'_manual-slab10.nii.gz',
                            output_dir=atlas_dir)['result']

    seg_result = nighres.registration.apply_coordinate_mappings(
                            image=seg_result,
                            mapping1=map_target, interpolation='nearest',
                            save_data=True, file_name=target+'_massp-slab10.nii.gz',
                            output_dir=atlas_dir)['result']

    nighres.statistics.segmentation_statistics(segmentation=seg_result, intensity=None, template=manual_img,
                                statistics=["Volumes","Volume_difference","Center_distance","Dice_overlap","Dilated_Dice_overlap","Average_surface_distance"], 
                                output_csv='atlas91_slab10_qmri2fcm_coregx2_overlap_v2_cl3_coreg2atlas.csv', file_name=target,
                                atlas=None, skip_first=True, ignore_zero=True)

