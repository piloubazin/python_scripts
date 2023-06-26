import nighres
import os
import nibabel
from glob import glob
import numpy
import scipy.ndimage

## use previously computed example as testbed
data_dir='/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/'
mapping_dir='/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/ahead-qmri2-mappings/'
in_dir ='/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/'
out_dir = in_dir+'qmri2_fcm/'
atlas_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/qmri2_ahead/'
reproc_dir ='/home/pilou/Projects/Ahead-Database/Automated-Parcellation/qmri2fcm_proc/'


subjects = ['sub-000','sub-001','sub-002','sub-004','sub-008','sub-012','sub-024','sub-031','sub-033','sub-040']
targets = ['sub-000','sub-001','sub-002','sub-004','sub-008','sub-012','sub-024','sub-031','sub-033','sub-040']
#targets = ['sub-001']

structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r',
              'ic_hem-l','ic_hem-r','vta_hem-l','vta_hem-r','fx_hem-lr','pag_hem-l','pag_hem-r','ppn_hem-l','ppn_hem-r',
              'cl_hem-l','cl_hem-r']


#target = targets[0]
for target in targets:

    n_contrasts = 4
    n_subjects = len(subjects)-1
    n_structures = len(structures)

    r1_target = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*r1hz_orient-std_brain.nii.gz')[0]
    r2s_target = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*r2hz_orient-std_brain.nii.gz')[0]
    qsm_target = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*qsm_orient-std_brain.nii.gz')[0]
    qpd_target = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*qpd_orient-std_brain.nii.gz')[0]
    
    # need to erode qsm mask
    qsm_img = nighres.io.load_volume(qsm_target)
    brainmask = (qsm_img.get_fdata()!=0)
    qsmmask = scipy.ndimage.binary_erosion(brainmask, iterations=5)
    qsm_img = nibabel.Nifti1Image(qsmmask*qsm_img.get_fdata(), qsm_img.affine, qsm_img.header)
    nighres.io.save_volume(atlas_dir+target+'_mod-qsm_masked.nii.gz', qsm_img)
    qsm_img = atlas_dir+target+'_mod-qsm_masked.nii.gz'
        
    # generate fcm-based contrasts
    rfcm = nighres.segmentation.fuzzy_cmeans(r1_target, clusters=2, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir)
    r1_target = rfcm['intensity']
    
    rfcm = nighres.segmentation.fuzzy_cmeans(r2s_target, clusters=2, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir)
    r2s_target = rfcm['intensity']
    
    rfcm = nighres.segmentation.fuzzy_cmeans(qsm_img, clusters=3, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir)
    qsm_target = rfcm['intensity']
    
    rfcm = nighres.segmentation.fuzzy_cmeans(qpd_target, clusters=3, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir)
    qpd_target = rfcm['intensity']

    map_target = glob(mapping_dir+target+'*map-ahead-qmri2_ants-map.nii.gz')[0]
    invmap_target = glob(mapping_dir+target+'*map-ahead-qmri2_ants-invmap.nii.gz')[0]
    
    manual_target = glob(in_dir+'label_maps/'+target+'_labeling-31.nii.gz')[0]
    
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
            r1_atlas = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2*r1hz_orient-std_brain.nii.gz')[0]
            r2s_atlas = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2*r2hz_orient-std_brain.nii.gz')[0]
            qsm_atlas = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2*qsm_orient-std_brain.nii.gz')[0]
            qpd_atlas = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2*qpd_orient-std_brain.nii.gz')[0]
            
            # need to erode qsm mask
            qsm_img = nighres.io.load_volume(qsm_atlas)
            brainmask = (qsm_img.get_fdata()!=0)
            qsmmask = scipy.ndimage.binary_erosion(brainmask, iterations=5)
            qsm_img = nibabel.Nifti1Image(qsmmask*qsm_img.get_fdata(), qsm_img.affine, qsm_img.header)
            nighres.io.save_volume(atlas_dir+sub+'_mod-qsm_masked.nii.gz', qsm_img)
            qsm_img = atlas_dir+sub+'_mod-qsm_masked.nii.gz'
            
            # generate fcm-based contrasts
            rfcm = nighres.segmentation.fuzzy_cmeans(r1_atlas, clusters=2, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir)
            r1_atlas = rfcm['intensity']
            
            rfcm = nighres.segmentation.fuzzy_cmeans(r2s_atlas, clusters=2, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir)
            r2s_atlas = rfcm['intensity']
            
            rfcm = nighres.segmentation.fuzzy_cmeans(qsm_img, clusters=3, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir)
            qsm_atlas = rfcm['intensity']
            
            rfcm = nighres.segmentation.fuzzy_cmeans(qpd_atlas, clusters=3, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir)
            qpd_atlas = rfcm['intensity']
            
            map_atlas = glob(mapping_dir+sub+'*map-ahead-qmri2_ants-map.nii.gz')[0]
            
            print("atlas subject: "+sub)
        
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=r1_atlas,
                            mapping1=map_atlas, interpolation='linear',
                            save_data=True, file_name=sub+'_mod-r1hz_def-qmri2fcm_atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
        
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=r2s_atlas,
                            mapping1=map_atlas, interpolation='linear',
                            save_data=True, file_name=sub+'_mod-r2hz_def-qmri2fcm_atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
        
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=qsm_atlas,
                            mapping1=map_atlas, interpolation='linear',
                            save_data=True, file_name=sub+'_mod-qsm_def-qmri2fcm_atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
            
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=qpd_atlas,
                            mapping1=map_atlas, interpolation='linear',
                            save_data=True, file_name=sub+'_mod-qpd_def-qmri2fcm_atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
            
            atlas_mappings[i] = map_atlas
        
            for j,struct in enumerate(structures):
                
                print('structure: '+struct+' (subject: '+sub+')')
                
                mask = glob(in_dir+'manual_delineations/best/'+sub+'*mask-'+struct+'*.nii.gz')[0]
                # flip the masks because of the change of canonical orientation
                mask_img = nighres.io.load_volume(mask)
                mask_img = nibabel.Nifti1Image(numpy.flip(mask_img.get_fdata(),axis=0),mask_img.affine,mask_img.header)

                defmask = nighres.registration.apply_coordinate_mappings(image=mask_img,
                                    mapping1=map_atlas, interpolation='linear',
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
                
            i = i+1;
    
    # this run is only to build the atlas
    nighres.parcellation.massp_atlasing(levelset_images=atlas_levelsets,
                            contrast_images=atlas_contrasts, 
                            skeleton_images=atlas_skeletons,
                            subjects=n_subjects, structures=n_structures,
                            contrasts=n_contrasts,
                            save_data=True, file_name='atlas91_31struct_qmri2fcm-'+target,
                            output_dir=atlas_dir, overwrite=False)
    
    # here we use the recomputed coregistration to the atlas (new version for the paper)
    invmap = reproc_dir+target+'_ses-1_reg2ahead-inverse_def-img.nii.gz'
    
    nighres.parcellation.massp(target_images=target_contrasts,
                            map_to_target=invmap, 
                            shape_atlas_probas=atlas_dir+'atlas91_31struct_qmri2fcm-'+target+'_massp-sproba.nii.gz', 
                            shape_atlas_labels=atlas_dir+'atlas91_31struct_qmri2fcm-'+target+'_massp-slabel.nii.gz',
                            intensity_atlas_hist=atlas_dir+'atlas91_31struct_qmri2fcm-'+target+'_massp-chist.nii.gz',
                            skeleton_atlas_probas=atlas_dir+'atlas91_31struct_qmri2fcm-'+target+'_massp-kproba.nii.gz', 
                            skeleton_atlas_labels=atlas_dir+'atlas91_31struct_qmri2fcm-'+target+'_massp-klabel.nii.gz',
                            max_iterations=120, max_difference=0.1,
                            intensity_prior=0.25,
                            save_data=True, file_name=target+'_atlas91_31struct_qmri2fcm_coregx2',
                            output_dir=out_dir, overwrite=False)
    
    seg_result = out_dir+target+'_atlas91_31struct_qmri2fcm_coregx2_massp-label.nii.gz'
    
    manual_img = nighres.io.load_volume(manual_target)
    manual_img = nibabel.Nifti1Image(numpy.flip(manual_img.get_fdata(),axis=0),manual_img.affine,manual_img.header)

    nighres.statistics.segmentation_statistics(segmentation=seg_result, intensity=None, template=manual_img,
                                statistics=["Volumes","Volume_difference","Center_distance","Dice_overlap","Dilated_Dice_overlap","Average_surface_distance"], 
                                output_csv='atlas91_31struct_qmri2fcm_coregx2_overlap.csv', file_name=target,
                                atlas=None, skip_first=True, ignore_zero=True)
