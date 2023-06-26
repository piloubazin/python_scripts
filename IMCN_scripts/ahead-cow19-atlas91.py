import nighres
import os
import nibabel
from glob import glob
import numpy
from PIL import Image

## use previously computed example as testbed
data_dir='/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/'
mapping_dir='/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/ahead-qmri2-mappings/'
arteries_dir='/home/pilou/Projects/Ahead-Database/Vascular-Modeling/process_qmri2/'
in_dir ='/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/'
out_dir = in_dir+'qmri2_ahead/'
atlas_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/qmri2_ahead/'

manart_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/manual_delineations/2DMIP/'
manart_lbls = ['ACA','PCA','ICA+MCA','AComA','PComA']

subjects = ['sub-000','sub-001','sub-002','sub-004','sub-008','sub-012','sub-024','sub-031','sub-033','sub-040']
targets = ['sub-000','sub-001','sub-002','sub-004','sub-008','sub-012','sub-024','sub-031','sub-033','sub-040']

structures = ['str_hem-l','str_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r',
              'ic_hem-l','ic_hem-r',
              'cl_hem-l','cl_hem-r','cow_hem-lr','inh_hem-lr']


#target = targets[0]
for target in targets:

    n_contrasts = 1
    n_subjects = len(subjects)-1
    n_structures = len(structures)

    r1_target = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2_mod-r1hz_orient-std_brain.nii.gz')[0]
#    r2s_target = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2_mod-r2hz_orient-std_brain.nii.gz')[0]
#    qsm_target = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2_mod-qsm_orient-std_brain.nii.gz')[0]
#    qpd_target = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2_mod-qpd_orient-std_brain.nii.gz')[0]
    art_target = glob(arteries_dir+target+'*mod-r1hz_orient-std_brain_mvf-pv.nii.gz')[0]
    
    map_target = glob(mapping_dir+target+'*map-ahead-qmri2_ants-map.nii.gz')[0]
    invmap_target = glob(mapping_dir+target+'*map-ahead-qmri2_ants-invmap.nii.gz')[0]
    
#    manual_target = glob(in_dir+'label_maps/'+target+'_labeling-32.nii.gz')[0]
    manual_target = in_dir+'label_maps/'+target+'_labeling-19.nii.gz'
    
#    target_contrasts = [r1_target,r2s_target,qsm_target,qpd_target,art_target]
#    target_contrasts = [r1_target,r2s_target,qsm_target,qpd_target]
    target_contrasts = [r1_target]
    
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
            r1_atlas = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2_mod-r1hz_orient-std_brain.nii.gz')[0]
#            r2s_atlas = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2_mod-r2hz_orient-std_brain.nii.gz')[0]
#            qsm_atlas = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2_mod-qsm_orient-std_brain.nii.gz')[0]
#            qpd_atlas = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2_mod-qpd_orient-std_brain.nii.gz')[0]
#            art_atlas = glob(arteries_dir+sub+'*mod-r1hz_orient-std_brain_mvf-pv.nii.gz')[0]
            
            map_atlas = glob(mapping_dir+sub+'*map-ahead-qmri2_ants-map.nii.gz')[0]
            
            print("atlas subject: "+sub)
        
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=r1_atlas,
                            mapping1=map_atlas, interpolation='linear',
                            save_data=True, file_name=sub+'_mod-r1hz_def-qmri2wb2_atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
        
#            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
#                            image=r2s_atlas,
#                            mapping1=map_atlas, interpolation='linear',
#                            save_data=True, file_name=sub+'_mod-r2hz_def-qmri2wb2_atlas.nii.gz',
#                            output_dir=atlas_dir)['result'])
        
#            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
#                            image=qsm_atlas,
#                            mapping1=map_atlas, interpolation='linear',
#                            save_data=True, file_name=sub+'_mod-qsm_def-qmri2wb2_atlas.nii.gz',
#                            output_dir=atlas_dir)['result'])
            
#            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
#                            image=qpd_atlas,
#                            mapping1=map_atlas, interpolation='linear',
#                            save_data=True, file_name=sub+'_mod-qpd_def-qmri2wb2_atlas.nii.gz',
#                            output_dir=atlas_dir)['result'])
                        
#            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
#                            image=art_atlas,
#                            mapping1=map_atlas, interpolation='linear',
#                            save_data=True, file_name=sub+'_mod-art_def-qmri2wb2_atlas.nii.gz',
#                            output_dir=atlas_dir)['result'])
            
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
    nighres.segmentation.conditional_shape_atlasing(levelset_images=atlas_levelsets,
                            contrast_images=atlas_contrasts, 
                            skeleton_images=atlas_skeletons,
                            subjects=n_subjects, structures=n_structures,
                            contrasts=n_contrasts,
                            save_data=True, file_name='atlas91_cow20-'+target,
                            output_dir=out_dir, overwrite=False)
    
    nighres.segmentation.conditional_shape(target_images=target_contrasts,
                            map_to_atlas=map_target, map_to_target=invmap_target, 
                            structures=n_structures, 
                            contrasts=n_contrasts, 
                            shape_atlas_probas=out_dir+'atlas91_cow20-'+target+'_cspmax-sproba.nii.gz', 
                            shape_atlas_labels=out_dir+'atlas91_cow20-'+target+'_cspmax-slabel.nii.gz',
                            intensity_atlas_hist=out_dir+'atlas91_cow20-'+target+'_cspmax-chist.nii.gz',
                            skeleton_atlas_probas=out_dir+'atlas91_cow20-'+target+'_cspmax-kproba.nii.gz', 
                            skeleton_atlas_labels=out_dir+'atlas91_cow20-'+target+'_cspmax-klabel.nii.gz',
                            max_iterations=120, max_difference=0.1,
                            save_data=True, file_name=target+'_atlas91_cow20',
                            output_dir=out_dir, overwrite=False)
    
    seg_result = out_dir+target+'_atlas91_cow20_cspmax-label.nii.gz'
    
#    manual_img = nighres.io.load_volume(manual_target)
#    manual_img = nibabel.Nifti1Image(numpy.flip(manual_img.get_fdata(),axis=0),manual_img.affine,manual_img.header)
#
#    nighres.statistics.segmentation_statistics(segmentation=seg_img, intensity=None, template=manual_img,
#                                statistics=["Volumes","Volume_difference","Center_distance","Dice_overlap","Dilated_Dice_overlap","Average_surface_distance"], 
#                                output_csv='atlas91_cow20_overlap.csv', file_name=target+'_stats',
#                                atlas=None, skip_first=True, ignore_zero=True)
    
    # multiply cow by vessel filter result
    art_data = nighres.io.load_volume(art_target).get_fdata()
    seg_img = nighres.io.load_volume(seg_result)
    
    seg_img = nibabel.Nifti1Image((seg_img.get_fdata()==19)*(art_data>0),seg_img.affine,seg_img.header)
    
    mask = glob(in_dir+'manual_delineations/best/'+target+'*mask-cow_hem-lr*.nii.gz')[0]
    mask_img = nighres.io.load_volume(mask)
    mask_img = nibabel.Nifti1Image(numpy.flip((mask_img.get_fdata()>0),axis=0),mask_img.affine,mask_img.header)                

#    nighres.statistics.segmentation_statistics(segmentation=seg_img, intensity=None, template=mask_img,
#                                statistics=["Volumes","Volume_difference","Center_distance","Dice_overlap","Dilated_Dice_overlap","Average_surface_distance"], 
#                                output_csv='atlas91_cow20xpv_overlap.csv', file_name=target+'_stats',
#                                atlas=None, skip_first=True, ignore_zero=True)
    
    # map manual to MNI
    mip_result = out_dir+target+'_cow-delineation_mip.png'
    if not os.path.exists(mip_result):
        defcow = nighres.registration.apply_coordinate_mappings(image=mask_img,
                mapping1=map_target, interpolation='linear',
                save_data=True, file_name=target+'_cow-manual', output_dir=out_dir)['result']

        # make a 2D MIP image
        defcow_data = nighres.io.load_volume(defcow).get_fdata()
        defcow_mip = numpy.transpose(numpy.max(defcow_data,axis=2))
        mip = Image.fromarray(numpy.uint8(defcow_mip*255))
        #mip = Image.fromarray(defcow_mip.astype('int'), mode='L')
        mip.save(mip_result)
            
    # map pv to MNI
#    mip_result = out_dir+target+'_cow-delineation_mip10.nii.gz'
#    if not os.path.exists(mip_result):
#        defcow = nighres.registration.apply_coordinate_mappings(image=mask_img,
#                mapping1=map_target, interpolation='linear',
#                save_data=True, file_name=target+'_cow-manual', output_dir=out_dir)['result']

        # make a partial MIP image
#        defcow_img = nighres.io.load_volume(defcow)
#        defcow_data = defcow_img.get_fdata()
        
#        mip = (defcow_data.shape[0],defcow_data.shape[1],defcow_data.shape[2]-10)
#        defcow_mip = numpy.zeros(mip)
#        for delta in range(10):
#            defcow_mip = numpy.maximum(defcow_mip,defcow_data[:,:,delta:-10+delta])
        
#        mip_img = nibabel.Nifti1Image(defcow_mip,defcow_img.affine,defcow_img.header)
#        nighres.io.save_volume(mip_result, mip_img)        
        
    # map manual arteries back
    for lbl in manart_lbls:
        manart = glob(manart_dir+target+'_2DMIP/*'+lbl+'*.png')
        manart_result = manart_dir+target+'_2DMIP/'+target+'_2dmip-'+lbl+'.nii.gz'
        if len(manart)>0 and not os.path.exists(manart_result):
            manart = manart[0]
            manart_data = numpy.transpose(numpy.asarray(Image.open(manart)))
            manart_data = manart_data.reshape(manart_data.shape[0],manart_data.shape[1],1)
            
            defcow = nighres.registration.apply_coordinate_mappings(image=mask_img,
                mapping1=map_target, interpolation='linear',
                save_data=True, file_name=target+'_cow-manual', output_dir=out_dir)['result']

            defcow_img = nighres.io.load_volume(defcow)
            defcow_data = defcow_img.get_fdata()*manart_data
            
            manart_img = nibabel.Nifti1Image(defcow_data,defcow_img.affine,defcow_img.header)
            nighres.io.save_volume(manart_result, manart_img)
            
            invcow = nighres.registration.apply_coordinate_mappings(image=manart_result,
                mapping1=invmap_target, interpolation='linear',
                save_data=True, output_dir=out_dir)['result']
            
            