import nighres
import os
import nibabel
from glob import glob
import numpy
import scipy.ndimage

## use previously computed example as testbed
data_dir='/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/'
in_dir ='/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/'
out_dir = '/home/pilou/Projects/Ahead-Database/Slab-Analysis/qmri2fcm_wb_slab_all/'
atlas_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/qmri2fcm_ahead/'
proc_dir ='/home/pilou/Projects/Ahead-Database/Automated-Parcellation/qmri2fcm_proc/'


wb_structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r',
              'ic_hem-l','ic_hem-r','vta_hem-l','vta_hem-r','fx_hem-lr','pag_hem-l','pag_hem-r','ppn_hem-l','ppn_hem-r',
              'cl_hem-l','cl_hem-r','ico_hem-l','ico_hem-r','sco_hem-l','sco_hem-r','lh_hem-l','lh_hem-r']

sb_structures = ['ac_hem-lr','pc_hem-lr','chn_hem-l','chn_hem-r','drn_hem-lr','mrn_hem-lr','rmg_hem-lr','rpo_hem-lr']

# regions that include each other
supregions = []
subregions = []   

session = 'ses-1'

for target_number in range(55,110): 
    target = 'sub-'+str(target_number).zfill(3)

    n_contrasts = 8
    n_structures = len(wb_structures)+len(sb_structures)

    r1_tsb = glob(data_dir+target+'/ses-1/anat/sb/qmri/'+target+'*sb2*r1hz_orient-std_brain.nii.gz')
    r2s_tsb = glob(data_dir+target+'/ses-1/anat/sb/qmri/'+target+'*sb2*r2hz_orient-std_brain.nii.gz')
    qsm_tsb = glob(data_dir+target+'/ses-1/anat/sb/qmri/'+target+'*sb2*qsm_orient-std_brain.nii.gz')
    qpd_tsb = glob(data_dir+target+'/ses-1/anat/sb/qmri/'+target+'*sb2*qpd_orient-std_brain.nii.gz')
    
    r1_twb = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*r1hz_orient-std_brain.nii.gz')
    r2s_twb = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*r2hz_orient-std_brain.nii.gz')
    qsm_twb = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*qsm_orient-std_brain.nii.gz')
    qpd_twb = glob(data_dir+target+'/ses-1/anat/wb/qmri/'+target+'*wb2*qpd_orient-std_brain.nii.gz')
    
    if len(r1_tsb)>0 and len(r2s_tsb)>0 and len(qsm_tsb)>0 and len(qpd_tsb)>0 and len(r1_twb)>0 and len(r2s_twb)>0 and len(qsm_twb)>0 and len(qpd_twb)>0:
        r1_tsb = r1_tsb[0]
        r2s_tsb = r2s_tsb[0]
        qsm_tsb = qsm_tsb[0]
        qpd_tsb = qpd_tsb[0]
        
        r1_twb = r1_twb[0]
        r2s_twb = r2s_twb[0]
        qsm_twb = qsm_twb[0]
        qpd_twb = qpd_twb[0]
    
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

        inv_tsb = nighres.registration.apply_coordinate_mappings(ants2['inverse'],mapping1=ants1['inverse'],mapping2=antsb1['mapping'],mapping3=antsb2['mapping'],
                                    save_data=True, overwrite=False, file_name=target+'_'+session+'_slab2ahead-inverse',
                                    output_dir=proc_dir)['result']
    
        inv_twb = nighres.registration.apply_coordinate_mappings(ants2['inverse'],mapping1=ants1['inverse'],
                                    save_data=True, overwrite=False, file_name=target+'_'+session+'_reg2ahead-inverse',
                                    output_dir=proc_dir)['result']
        
        del antsb1,antsb2,ants1,ants2
        del mapping_swb,mapping_twb
        

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
    
        parcell = nighres.parcellation.massp(target_images=target_contrasts,
                                map_to_target=None, structures=n_structures,
                                shape_atlas_probas=atlas_dir+'atlas10_wb-slab-raphe_nonzero_qmri2fcm_cl3_coregx4clean_hist5_cspmax-sproba.nii.gz', 
                                shape_atlas_labels=atlas_dir+'atlas10_wb-slab-raphe_nonzero_qmri2fcm_cl3_coregx4clean_hist5_cspmax-slabel.nii.gz',
                                intensity_atlas_hist=atlas_dir+'atlas10_wb-slab-raphe_nonzero_qmri2fcm_cl3_coregx4clean_hist5_cspmax-chist.nii.gz',
                                skeleton_atlas_probas=atlas_dir+'atlas10_wb-slab-raphe_nonzero_qmri2fcm_cl3_coregx4clean_hist5_cspmax-kproba.nii.gz', 
                                skeleton_atlas_labels=atlas_dir+'atlas10_wb-slab-raphe_nonzero_qmri2fcm_cl3_coregx4clean_hist5_cspmax-klabel.nii.gz',
                                max_iterations=120, max_difference=0.1,
                                intensity_prior=0.25,
                                save_data=True, file_name=target+'_atlas10_wb-slab-raphe_nonzero_qmri2fcm_cl3_coregx4clean_hist5',
                                output_dir=out_dir, overwrite=False)
    
        # map the parcellations back into whole brain and slab space
        nighres.registration.apply_coordinate_mappings(
                                image=parcell['max_label'],
                                mapping1=inv_tsb, interpolation='nearest', padding='zero', zero_border=0,
                                save_data=True, overwrite=False, file_name=target+'_atlas10_wb-slab-raphe_nonzero_qmri2fcm_cl3_coregx4clean_hist5_massp-label_to-sb.nii.gz',
                                output_dir=out_dir)
        
        nighres.registration.apply_coordinate_mappings(
                                image=parcell['max_label'],
                                mapping1=inv_twb, interpolation='nearest', padding='zero', zero_border=0,
                                save_data=True, overwrite=False, file_name=target+'_atlas10_wb-slab-raphe_nonzero_qmri2fcm_cl3_coregx4clean_hist5_massp-label_to-wb.nii.gz',
                                output_dir=out_dir)
        
        
        
