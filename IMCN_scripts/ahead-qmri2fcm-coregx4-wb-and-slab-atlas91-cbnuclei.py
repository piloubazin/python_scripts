import nighres
import os
import nibabel
from glob import glob
import numpy
import scipy.ndimage

## use previously computed example as testbed
data_dir='/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/'
in_dir ='/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/'
out_dir = '/home/pilou/Projects/Ahead-Database/Slab-Analysis/qmri2fcm_wb_slab/'
atlas_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/qmri2fcm_ahead/'
proc_dir ='/home/pilou/Projects/Ahead-Database/Automated-Parcellation/qmri2fcm_proc/'


subjects = ['sub-000','sub-001','sub-002','sub-004','sub-008','sub-012','sub-024','sub-031','sub-033','sub-040']
targets = ['sub-000','sub-001','sub-002','sub-004','sub-008','sub-012','sub-024','sub-031','sub-033','sub-040']
#targets = ['sub-033']

wb_structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r',
              'ic_hem-l','ic_hem-r','vta_hem-l','vta_hem-r','fx_hem-lr','pag_hem-l','pag_hem-r','ppn_hem-l','ppn_hem-r',
              'cl_hem-l','cl_hem-r','ico_hem-l','ico_hem-r','sco_hem-l','sco_hem-r','lh_hem-l','lh_hem-r']

sb_structures = ['ac_hem-lr','pc_hem-lr','chn_hem-l','chn_hem-r','drn_hem-lr',
                'cbdn_hem-l','cbdn_hem-r','cbem_hem-l','cbem_hem-r','cbgl_hem-l','cbgl_hem-r','cbfa_hem-lr']

# regions that include each other
supregions = []
subregions = []   

session = 'ses-1'

for target in targets:

    n_contrasts = 8
    n_subjects = len(subjects)-1
    n_structures = len(wb_structures)+len(sb_structures)

    r1_tsb = glob(data_dir+target+'/ses-1/anat/sb/qmri/'+target+'*sb2*r1hz_orient-std_brain.nii.gz')[0]
    r2s_tsb = glob(data_dir+target+'/ses-1/anat/sb/qmri/'+target+'*sb2*r2hz_orient-std_brain.nii.gz')[0]
    qsm_tsb = glob(data_dir+target+'/ses-1/anat/sb/qmri/'+target+'*sb2*qsm_orient-std_brain.nii.gz')[0]
    qpd_tsb = glob(data_dir+target+'/ses-1/anat/sb/qmri/'+target+'*sb2*qpd_orient-std_brain.nii.gz')[0]
    
    r1_twb = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*r1hz_orient-std_brain.nii.gz')[0]
    r2s_twb = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*r2hz_orient-std_brain.nii.gz')[0]
    qsm_twb = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*qsm_orient-std_brain.nii.gz')[0]
    qpd_twb = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*qpd_orient-std_brain.nii.gz')[0]
    
    # generate fcm-based contrasts: maybe not a good idea here?
    rfcm = nighres.segmentation.fuzzy_cmeans(r1_tsb, clusters=3, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir, 
                                             file_name=target+'_mod-r1hz_clusters-3.nii.gz')
    r1_tsb = rfcm['intensity']
    
    rfcm = nighres.segmentation.fuzzy_cmeans(r2s_tsb, clusters=3, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir, 
                                             file_name=target+'_mod-r2hz_clusters-3.nii.gz')
    r2s_tsb = rfcm['intensity']
    
    rfcm = nighres.segmentation.fuzzy_cmeans(qsm_tsb, clusters=3, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir)
    qsm_tsb = rfcm['intensity']
    
    rfcm = nighres.segmentation.fuzzy_cmeans(qpd_tsb, clusters=3, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir)
    qpd_tsb = rfcm['intensity']


    rfcm = nighres.segmentation.fuzzy_cmeans(r1_twb, clusters=3, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir, 
                                             file_name=target+'_wb_mod-r1hz_clusters-3.nii.gz')
    r1_twb = rfcm['intensity']
    
    rfcm = nighres.segmentation.fuzzy_cmeans(r2s_twb, clusters=3, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir, 
                                             file_name=target+'_wb_mod-r2hz_clusters-3.nii.gz')
    r2s_twb = rfcm['intensity']
    
    rfcm = nighres.segmentation.fuzzy_cmeans(qsm_twb, clusters=3, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir)
    qsm_twb = rfcm['intensity']
    
    rfcm = nighres.segmentation.fuzzy_cmeans(qpd_twb, clusters=3, max_iterations=150, max_difference=0.001, 
                                             smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                             save_data=True, overwrite=False, output_dir=atlas_dir)
    qpd_twb = rfcm['intensity']

    del rfcm

    # explicitly recompute a better coreg
    template = nighres.data.download_AHEAD_template()
    
    # first computing the coreg to whole brain
    antsb1 = nighres.registration.embedded_antspy_multi(
                                    source_images=[r1_twb,r2s_twb],
                                    target_images=[r1_tsb,r2s_tsb],
                                    run_rigid=True, run_affine=False, run_syn=False,
                                    rigid_iterations=10000,
                                    affine_iterations=2000,
                                    coarse_iterations=180, 
                                    medium_iterations=60, fine_iterations=30,
                                    cost_function='MutualInformation', 
                                    interpolation='Linear',
                                    regularization='High',
                                    smooth_mask=0.1,
                                    ignore_affine=False, 
                                    save_data=True, file_name=target+'_'+session+'_reg2slab-step1',
                                    output_dir=proc_dir)

    antsb2 = nighres.registration.embedded_antspy_multi(
                                    source_images=[antsb1['transformed_sources'][0],antsb1['transformed_sources'][1]],
                                    target_images=[r1_tsb,r2s_tsb],
                                    run_rigid=True, run_affine=False, run_syn=True,
                                    rigid_iterations=10000,
                                    affine_iterations=2000,
                                    coarse_iterations=180, 
                                    medium_iterations=60, fine_iterations=30,
                                    cost_function='MutualInformation', 
                                    interpolation='Linear',
                                    regularization='High',
                                    smooth_mask=0.1,
                                    ignore_affine=False, 
                                    save_data=True, file_name=target+'_'+session+'_reg2slab-step2',
                                    output_dir=proc_dir)

    ants1 = nighres.registration.embedded_antspy_multi(
                                    source_images=[r1_twb,r2s_twb],
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

    mapping_swb = nighres.registration.apply_coordinate_mappings(antsb2['inverse'],mapping1=antsb1['inverse'],mapping2=ants1['mapping'],mapping3=ants2['mapping'],
                                save_data=True, overwrite=False, file_name=target+'_'+session+'_slab2ahead-combined',
                                output_dir=proc_dir)

    mapping_twb = nighres.registration.apply_coordinate_mappings(ants1['mapping'],mapping1=ants2['mapping'],
                                save_data=True, overwrite=False, file_name=target+'_'+session+'_reg2ahead-mapping',
                                output_dir=proc_dir)
    
    map_tsb = mapping_swb['result']
    map_twb = mapping_twb['result']
    
    del antsb1,antsb2,ants1,ants2
    del mapping_swb,mapping_twb
    
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
            r1_asb = glob(data_dir+sub+'/ses-1/anat/sb/qmri/'+sub+'*sb2*r1hz_orient-std_brain.nii.gz')[0]
            r2s_asb = glob(data_dir+sub+'/ses-1/anat/sb/qmri/'+sub+'*sb2*r2hz_orient-std_brain.nii.gz')[0]
            qsm_asb = glob(data_dir+sub+'/ses-1/anat/sb/qmri/'+sub+'*sb2*qsm_orient-std_brain.nii.gz')[0]
            qpd_asb = glob(data_dir+sub+'/ses-1/anat/sb/qmri/'+sub+'*sb2*qpd_orient-std_brain.nii.gz')[0]
            
            r1_awb = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2*r1hz_orient-std_brain.nii.gz')[0]
            r2s_awb = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2*r2hz_orient-std_brain.nii.gz')[0]
            qsm_awb = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2*qsm_orient-std_brain.nii.gz')[0]
            qpd_awb = glob(data_dir+sub+'/ses-1/anat/wb/qmri/'+sub+'*wb2*qpd_orient-std_brain.nii.gz')[0]
            
            # generate fcm-based contrasts
            rfcm = nighres.segmentation.fuzzy_cmeans(r1_asb, clusters=3, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir, 
                                                     file_name=sub+'_mod-r1hz_clusters-3.nii.gz')
            r1_asb = rfcm['intensity']
            
            rfcm = nighres.segmentation.fuzzy_cmeans(r2s_asb, clusters=2, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir, 
                                                     file_name=sub+'_mod-r2hz_clusters-3.nii.gz')
            r2s_asb = rfcm['intensity']
            
            rfcm = nighres.segmentation.fuzzy_cmeans(qsm_asb, clusters=3, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir)
            qsm_asb = rfcm['intensity']
            
            rfcm = nighres.segmentation.fuzzy_cmeans(qpd_asb, clusters=3, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir)
            qpd_asb = rfcm['intensity']


            rfcm = nighres.segmentation.fuzzy_cmeans(r1_awb, clusters=3, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir, 
                                                     file_name=sub+'_wb_mod-r1hz_clusters-3.nii.gz')
            r1_awb = rfcm['intensity']
            
            rfcm = nighres.segmentation.fuzzy_cmeans(r2s_awb, clusters=3, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir, 
                                                     file_name=sub+'_wb_mod-r2hz_clusters-3.nii.gz')
            r2s_awb = rfcm['intensity']
            
            rfcm = nighres.segmentation.fuzzy_cmeans(qsm_awb, clusters=3, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir)
            qsm_awb = rfcm['intensity']
            
            rfcm = nighres.segmentation.fuzzy_cmeans(qpd_awb, clusters=3, max_iterations=150, max_difference=0.001, 
                                                     smoothing=0.0, fuzziness=2.0, mask_zero=True, map_intensity=True,
                                                     save_data=True, overwrite=False, output_dir=atlas_dir)
            qpd_awb = rfcm['intensity']

            del rfcm

            antsb1 = nighres.registration.embedded_antspy_multi(
                                            source_images=[r1_awb,r2s_awb],
                                            target_images=[r1_asb,r2s_asb],
                                            run_rigid=True, run_affine=False, run_syn=False,
                                            rigid_iterations=10000,
                                            affine_iterations=2000,
                                            coarse_iterations=180, 
                                            medium_iterations=60, fine_iterations=30,
                                            cost_function='MutualInformation', 
                                            interpolation='Linear',
                                            regularization='High',
                                            smooth_mask=0.1,
                                            ignore_affine=False, 
                                            save_data=True, file_name=sub+'_'+session+'_reg2slab-step1',
                                            output_dir=proc_dir)
        
            antsb2 = nighres.registration.embedded_antspy_multi(
                                            source_images=[antsb1['transformed_sources'][0],antsb1['transformed_sources'][1]],
                                            target_images=[r1_asb,r2s_asb],
                                            run_rigid=True, run_affine=False, run_syn=True,
                                            rigid_iterations=10000,
                                            affine_iterations=2000,
                                            coarse_iterations=180, 
                                            medium_iterations=60, fine_iterations=30,
                                            cost_function='MutualInformation', 
                                            interpolation='Linear',
                                            regularization='High',
                                            smooth_mask=0.1,
                                            ignore_affine=False, 
                                            save_data=True, file_name=sub+'_'+session+'_reg2slab-step2',
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
        
            mapping_asb = nighres.registration.apply_coordinate_mappings(antsb2['inverse'],mapping1=antsb1['inverse'],mapping2=ants1['mapping'],mapping3=ants2['mapping'],
                                        save_data=True, overwrite=False, file_name=sub+'_'+session+'_slab2ahead-combined',
                                        output_dir=proc_dir)

            mapping_awb = nighres.registration.apply_coordinate_mappings(ants1['mapping'],mapping1=ants2['mapping'],
                                        save_data=True, overwrite=False, file_name=sub+'_'+session+'_reg2ahead-mapping',
                                        output_dir=proc_dir)
        
            map_asb = mapping_asb['result']
            map_awb = mapping_awb['result']
            
            del ants1,ants2,antsb1,antsb2
            del mapping_asb,mapping_awb
            
            print("atlas subject: "+sub)
        
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=r1_awb,
                            mapping1=map_awb, interpolation='linear',
                            save_data=True, overwrite=False, file_name=sub+'_mod-r1hz_def-qmri2fcm_cl3-atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
        
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=r2s_awb,
                            mapping1=map_awb, interpolation='linear',
                            save_data=True, overwrite=False, file_name=sub+'_mod-r2hz_def-qmri2fcm_cl3-atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
        
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=qsm_awb,
                            mapping1=map_awb, interpolation='linear',
                            save_data=True, file_name=sub+'_mod-qsm_def-qmri2fcm_atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
            
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=qpd_awb,
                            mapping1=map_awb, interpolation='linear',
                            save_data=True, file_name=sub+'_mod-qpd_def-qmri2fcm_atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
            
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=r1_asb,
                            mapping1=map_asb, interpolation='linear', padding='zero', zero_border=10,
                            save_data=True, overwrite=False, file_name=sub+'_mod-r1hz_defx4-qmri2fcm_cl3-slab-atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
        
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=r2s_asb,
                            mapping1=map_asb, interpolation='linear', padding='zero', zero_border=10,
                            save_data=True, overwrite=False, file_name=sub+'_mod-r2hz_defx4-qmri2fcm_cl3-slab-atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
        
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=qsm_asb,
                            mapping1=map_asb, interpolation='linear', padding='zero', zero_border=10,
                            save_data=True, overwrite=False, file_name=sub+'_mod-qsm_defx4-qmri2fcm_slab-atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
            
            atlas_contrasts[i].append(nighres.registration.apply_coordinate_mappings(
                            image=qpd_asb,
                            mapping1=map_asb, interpolation='linear', padding='zero', zero_border=10,
                            save_data=True, overwrite=False, file_name=sub+'_mod-qpd_defx4-qmri2fcm_slab-atlas.nii.gz',
                            output_dir=atlas_dir)['result'])
            
            atlas_mappings[i] = map_asb
        
            for j,struct in enumerate(wb_structures):
                
                print('structure: '+struct+' (subject: '+sub+')')
                
                mask = glob(in_dir+'manual_delineations/best/'+sub+'*_mask-'+struct+'*.nii.gz')[0]
                # flip the masks because of the change of canonical orientation
                # only to be done for old masks, check sources for new ones, if any
                mask_img = nighres.io.load_volume(mask)
                mask_img = nibabel.Nifti1Image(numpy.flip(mask_img.get_fdata(),axis=0),mask_img.affine,mask_img.header)

                defmask = nighres.io.load_volume(nighres.registration.apply_coordinate_mappings(image=mask_img,
                                    mapping1=map_awb, interpolation='linear',
                                    save_data=True, file_name=sub+'_mask-'+struct+'_def-qmri2_atlas.nii.gz',
                                    output_dir=atlas_dir)['result'])
                
                if struct in supregions:
                    for idx,supreg in enumerate(supregions):
                        if supreg==struct:
                            
                            subreg=subregions[idx]
                            structname = struct
                            
                            roi = nighres.io.load_volume(glob(in_dir+'manual_delineations/best/'+sub+'*_mask-'+subreg+'_*.nii.gz')[0])
                            roi = nibabel.as_closest_canonical(roi)
                            roi = nighres.io.load_volume(nighres.registration.apply_coordinate_mappings(image=roi,
                                    mapping1=map_asb, interpolation='linear',
                                    save_data=True, file_name=sub+'_mask-'+subreg+'_def-qmri2_slab-atlas.nii.gz',
                                    output_dir=atlas_dir)['result']).get_fdata()
                            
                            roi = numpy.maximum(0.0, defmask.get_fdata()-roi)
                
                            defmask = nibabel.Nifti1Image(roi,mask_img.affine,mask_img.header)
                            
                            structname = structname+'_sub-'+subreg
                            
                            del roi
                            
                    struct = structname            
                
                
                atlas_levelsets[i].append(nighres.surface.probability_to_levelset(
                                        probability_image=defmask,
                                        save_data=True, file_name=sub+'_mask-'+struct+'_def-qmri2_atlas.nii.gz',
                                        output_dir=atlas_dir)['result'])
                atlas_skeletons[i].append(nighres.shape.levelset_thickness(
                                        atlas_levelsets[i][j],
                                        save_data=True, file_name=sub+'_mask-'+struct+'_def-qmri2_atlas.nii.gz',
                                        output_dir=atlas_dir)['dist'])
                
                del mask_img,defmask
                
            for j,struct in enumerate(sb_structures):
                
                print('structure: '+struct+' (subject: '+sub+')')
                
                # check if exists: otherwise look for the whole brain version
                mask = glob(in_dir+'manual_delineations/best/'+sub+'_ses-1_acq-sb*_mask-'+struct+'*.nii.gz')
                if len(mask)>0:
                
                    mask = glob(in_dir+'manual_delineations/best/'+sub+'_ses-1_acq-sb*_mask-'+struct+'*.nii.gz')[0]
                    # do NOT flip the masks that were done on the new slab data!!
                    mask_img = nighres.io.load_volume(mask)
                    #mask_img = nibabel.Nifti1Image(numpy.flip(mask_img.get_fdata(),axis=0),mask_img.affine,mask_img.header)
                    
                    # note: some slabs were oriented incorrectly in voxel space: align first
                    mask_img = nibabel.as_closest_canonical(mask_img)
    
                    defmask = nighres.registration.apply_coordinate_mappings(image=mask_img,
                                        mapping1=map_asb, interpolation='linear',
                                        save_data=True, overwrite=False, file_name=sub+'_mask-'+struct+'_def-qmri2_slab-atlas.nii.gz',
                                        output_dir=atlas_dir)['result']
                    atlas_levelsets[i].append(nighres.surface.probability_to_levelset(
                                            probability_image=defmask,
                                            save_data=True, overwrite=False, file_name=sub+'_mask-'+struct+'_def-qmri2_slab-atlas.nii.gz',
                                            output_dir=atlas_dir)['result'])
                    atlas_skeletons[i].append(nighres.shape.levelset_thickness(
                                            atlas_levelsets[i][j],
                                            save_data=True, overwrite=False, file_name=sub+'_mask-'+struct+'_def-qmri2_slab-atlas.nii.gz',
                                            output_dir=atlas_dir)['dist'])
                    
                    del mask_img,defmask
                    
                else:
                    mask = glob(in_dir+'manual_delineations/best/'+sub+'_ses-1_acq-wb*_mask-'+struct+'*.nii.gz')[0]
                    
                    # do NOT flip the masks that were done on the new slab data!!
                    mask_img = nighres.io.load_volume(mask)
                    #mask_img = nibabel.Nifti1Image(numpy.flip(mask_img.get_fdata(),axis=0),mask_img.affine,mask_img.header)
                    
                    # note: some slabs were oriented incorrectly in voxel space: align first
                    mask_img = nibabel.as_closest_canonical(mask_img)
    
                    defmask = nighres.registration.apply_coordinate_mappings(image=mask_img,
                                        mapping1=map_awb, interpolation='linear',
                                        save_data=True, overwrite=False, file_name=sub+'_mask-'+struct+'_def-qmri2_slab-atlas.nii.gz',
                                        output_dir=atlas_dir)['result']
                    atlas_levelsets[i].append(nighres.surface.probability_to_levelset(
                                            probability_image=defmask,
                                            save_data=True, overwrite=False, file_name=sub+'_mask-'+struct+'_def-qmri2_slab-atlas.nii.gz',
                                            output_dir=atlas_dir)['result'])
                    atlas_skeletons[i].append(nighres.shape.levelset_thickness(
                                            atlas_levelsets[i][j],
                                            save_data=True, overwrite=False, file_name=sub+'_mask-'+struct+'_def-qmri2_slab-atlas.nii.gz',
                                            output_dir=atlas_dir)['dist'])
                    
                    del mask_img,defmask
                    
                
            i = i+1;
    
    # this run is only to build the atlas
    nighres.segmentation.conditional_shape_atlasing(levelset_images=atlas_levelsets,
                            contrast_images=atlas_contrasts, 
                            skeleton_images=atlas_skeletons,
                            subjects=n_subjects, structures=n_structures,
                            contrasts=n_contrasts, smoothing=5.0,
                            save_data=True, file_name='atlas91_wb-slab-cb_nonzero_qmri2fcm_cl3_coregx4clean_hist5-'+target,
                            output_dir=atlas_dir, overwrite=False)


    # bring all subject data in MNI space
    r1_twb = nighres.registration.apply_coordinate_mappings(
                            image=r1_twb,
                            mapping1=map_twb, interpolation='linear',
                            save_data=True, overwrite=False, file_name=target+'_mod-r1hz_def-qmri2fcm_cl3-atlas.nii.gz',
                            output_dir=atlas_dir)['result']
        
    r2s_twb = nighres.registration.apply_coordinate_mappings(
                            image=r2s_twb,
                            mapping1=map_twb, interpolation='linear',
                            save_data=True, overwrite=False, file_name=target+'_mod-r2hz_def-qmri2fcm_cl3-atlas.nii.gz',
                            output_dir=atlas_dir)['result']
        
    qsm_twb = nighres.registration.apply_coordinate_mappings(
                            image=qsm_twb,
                            mapping1=map_twb, interpolation='linear',
                            save_data=True, file_name=target+'_mod-qsm_def-qmri2fcm_atlas.nii.gz',
                            output_dir=atlas_dir)['result']
            
    qpd_twb = nighres.registration.apply_coordinate_mappings(
                            image=qpd_twb,
                            mapping1=map_twb, interpolation='linear',
                            save_data=True, file_name=target+'_mod-qpd_def-qmri2fcm_atlas.nii.gz',
                            output_dir=atlas_dir)['result']
            
    r1_tsb = nighres.registration.apply_coordinate_mappings(
                            image=r1_tsb,
                            mapping1=map_tsb, interpolation='linear', padding='zero', zero_border=10,
                            save_data=True, overwrite=False, file_name=target+'_mod-r1hz_defx4-qmri2fcm_cl3-slab-atlas.nii.gz',
                            output_dir=atlas_dir)['result']
        
    r2s_tsb = nighres.registration.apply_coordinate_mappings(
                            image=r2s_tsb,
                            mapping1=map_tsb, interpolation='linear', padding='zero', zero_border=10,
                            save_data=True, overwrite=False, file_name=target+'_mod-r2hz_defx4-qmri2fcm_cl3-slab-atlas.nii.gz',
                            output_dir=atlas_dir)['result']
        
    qsm_tsb = nighres.registration.apply_coordinate_mappings(
                            image=qsm_tsb,
                            mapping1=map_tsb, interpolation='linear', padding='zero', zero_border=10,
                            save_data=True, overwrite=False, file_name=target+'_mod-qsm_defx4-qmri2fcm_slab-atlas.nii.gz',
                            output_dir=atlas_dir)['result']
            
    qpd_tsb = nighres.registration.apply_coordinate_mappings(
                            image=qpd_tsb,
                            mapping1=map_tsb, interpolation='linear', padding='zero', zero_border=10,
                            save_data=True, overwrite=False, file_name=target+'_mod-qpd_defx4-qmri2fcm_slab-atlas.nii.gz',
                            output_dir=atlas_dir)['result']

    target_contrasts = [r1_twb,r2s_twb,qsm_twb,qpd_twb,r1_tsb,r2s_tsb,qsm_tsb,qpd_tsb]

    nighres.parcellation.massp(target_images=target_contrasts,
                            map_to_target=None, structures=n_structures,
                            shape_atlas_probas=atlas_dir+'atlas91_wb-slab-cb_nonzero_qmri2fcm_cl3_coregx4clean_hist5-'+target+'_cspmax-sproba.nii.gz', 
                            shape_atlas_labels=atlas_dir+'atlas91_wb-slab-cb_nonzero_qmri2fcm_cl3_coregx4clean_hist5-'+target+'_cspmax-slabel.nii.gz',
                            intensity_atlas_hist=atlas_dir+'atlas91_wb-slab-cb_nonzero_qmri2fcm_cl3_coregx4clean_hist5-'+target+'_cspmax-chist.nii.gz',
                            skeleton_atlas_probas=atlas_dir+'atlas91_wb-slab-cb_nonzero_qmri2fcm_cl3_coregx4clean_hist5-'+target+'_cspmax-kproba.nii.gz', 
                            skeleton_atlas_labels=atlas_dir+'atlas91_wb-slab-cb_nonzero_qmri2fcm_cl3_coregx4clean_hist5-'+target+'_cspmax-klabel.nii.gz',
                            max_iterations=120, max_difference=0.1,
                            intensity_prior=0.25,
                            save_data=True, file_name=target+'_atlas91_wb-slab-cb_nonzero_qmri2fcm_cl3_coregx4clean_hist5',
                            output_dir=out_dir, overwrite=False)
    
    seg_result = out_dir+target+'_atlas91_wb-slab-cb_nonzero_qmri2fcm_cl3_coregx4clean_hist5_massp-label.nii.gz'
    
    manual_target = glob(in_dir+'label_maps/'+target+'_labeling-wb-and-slab49.nii.gz')[0]
    
    #manual_img = nighres.io.load_volume(manual_target)
    #manual_img = nibabel.Nifti1Image(numpy.flip(manual_img.get_fdata(),axis=0),manual_img.affine,manual_img.header)

    nighres.statistics.segmentation_statistics(segmentation=seg_result, intensity=None, template=manual_target,
                                statistics=["Volumes","Volume_difference","Center_distance","Dice_overlap","Dilated_Dice_overlap","Average_surface_distance"], 
                                output_csv='atlas91_wb-slab-cb_nonzero_qmri2fcm_cl3_coregx4clean_hist5_overlap.csv', file_name=target,
                                atlas=None, skip_first=True, ignore_zero=True)

