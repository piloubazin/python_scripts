import nibabel
import nighres
from glob import glob
import numpy
import os
from nighres.io import load_volume, save_volume

delineation_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/manual_delineations/best'
labeling_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/label_maps'
mesh_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/meshes'
data_dir='/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/'
atlas_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/qmri2fcm_ahead/'
proc_dir ='/home/pilou/Projects/Ahead-Database/Automated-Parcellation/qmri2fcm_proc/'


subjects = ['sub-000','sub-001','sub-002','sub-004','sub-008','sub-012','sub-024','sub-031','sub-033','sub-040']
#subjects = ['sub-033']

wb_structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r',
              'ic_hem-l','ic_hem-r','vta_hem-l','vta_hem-r','fx_hem-lr','pag_hem-l','pag_hem-r','ppn_hem-l','ppn_hem-r',
              'cl_hem-l','cl_hem-r','ico_hem-l','ico_hem-r','sco_hem-l','sco_hem-r','lh_hem-l','lh_hem-r']

#sb_structures = ['ac_hem-lr','pc_hem-lr','chn_hem-l','chn_hem-r','drn_hem-lr',
#              'lgn_hem-l','lgn_hem-r','mgn_hem-l','mgn_hem-r','ns_hem-l','ns_hem-r']
#                'cbdn_hem-l','cbdn_hem-r','cbem_hem-l','cbem_hem-r','cbgl_hem-l','cbgl_hem-r','cbfa_hem-lr']
sb_structures = ['ac_hem-lr','pc_hem-lr','chn_hem-l','chn_hem-r','drn_hem-lr','mrn_hem-lr','rmg_hem-lr','rpo_hem-lr']

# regions that include each other
#supregions = ['tha_hem-l','tha_hem-r','tha_hem-l','tha_hem-r','sn_hem-l','sn_hem-r']
#subregions = ['lgn_hem-l','lgn_hem-r','mgn_hem-l','mgn_hem-r','ns_hem-l','ns_hem-r']   
supregions = []
subregions = []


session = 'ses-1'


build_mesh = False

label_map = None
min_dist = None
for subject in subjects:
    # find transformation files
    ants0 = proc_dir+subject+'_'+session+'_reg2slab-step0_ants-map.nii.gz'

    ants1 = proc_dir+subject+'_'+session+'_reg2ahead-step1_ants-map.nii.gz'
            
    ants2 = proc_dir+subject+'_'+session+'_reg2ahead-step2_ants-map.nii.gz'

    map_sb = nighres.registration.apply_coordinate_mappings(ants0,mapping1=ants1,mapping2=ants2,
                                save_data=True, overwrite=False, file_name=subject+'_'+session+'_slab2ahead-mapping',
                                output_dir=proc_dir)['result']

    map_wb = nighres.registration.apply_coordinate_mappings(ants1,mapping1=ants2,
                                save_data=True, overwrite=False, file_name=subject+'_'+session+'_reg2ahead-mapping',
                                output_dir=proc_dir)['result']
    
    label_file = labeling_dir+'/'+subject+'_labeling-wb-and-slab'+str(len(wb_structures)+len(sb_structures))+'.nii.gz'
    
    if not os.path.exists(label_file):
        label = 1
        for structure in wb_structures:
                
            # find the first instance of segmentation
            segfiles = glob(delineation_dir+'/'+subject+'*_mask-'+structure+'_*.nii.gz')
            
            seg = None
            if len(segfiles)>0:
                seg = load_volume(segfiles[0])
                
            if seg is not None:    
                print('load :'+segfiles[0])
                
                #seg = nibabel.as_closest_canonical(seg)
                seg = nibabel.Nifti1Image(numpy.flip(seg.get_fdata(),axis=0),seg.affine,seg.header)
                #nighres.io.save_volume(labeling_dir+'/'+subject+'_reorient_mask-'+structure+'.nii.gz',seg)
                
                seg = nighres.io.load_volume(nighres.registration.apply_coordinate_mappings(image=seg,
                                    mapping1=map_wb, interpolation='linear',
                                    save_data=True, file_name=subject+'_mask-'+structure+'_def-qmri2_atlas.nii.gz',
                                    output_dir=atlas_dir)['result'])
                
                if structure in supregions:
                    for idx,supreg in enumerate(supregions):
                        if supreg==structure:
                            
                            subreg=subregions[idx]
                            structname = structure
                            
                            roi = load_volume(glob(delineation_dir+'/'+subject+'*_mask-'+subreg+'_*.nii.gz')[0])
                            roi = nibabel.as_closest_canonical(roi)
                            roi = nighres.io.load_volume(nighres.registration.apply_coordinate_mappings(image=roi,
                                    mapping1=map_sb, interpolation='linear',
                                    save_data=True, file_name=subject+'_mask-'+subreg+'_def-qmri2_slab-atlas.nii.gz',
                                    output_dir=atlas_dir)['result']).get_fdata()
                            
                            roi = numpy.maximum(0.0, seg.get_fdata()-roi)
                
                            seg = nibabel.Nifti1Image(roi,seg.affine,seg.header)
                            
                            structname = structname+'_sub-'+subreg
                            
                    structure = structname            

                lvlset = nighres.io.load_volume(nighres.surface.probability_to_levelset(probability_image=seg,
                                        save_data=True, file_name=subject+'_mask-'+structure+'_def-qmri2_atlas.nii.gz',
                                        output_dir=atlas_dir)['result']).get_fdata()
                
                if label_map is None:
                    min_dist = lvlset
                    label_map = nibabel.Nifti1Image((lvlset<=0)*label, seg.affine, seg.header)
                    min_dist[(lvlset>0)] = 0
                else:
                    data = numpy.maximum(label_map.get_fdata(),(lvlset<min_dist)*label)
                    label_map = nibabel.Nifti1Image(data, seg.affine, seg.header)
                    min_dist = numpy.minimum(min_dist,lvlset)
                
                if build_mesh:
                    mesh_file = subject+'_mesh-'+structure+'.vtk'
                    mesh = nighres.surface.levelset_to_mesh(lvlset,save_data=True,output_dir=mesh_dir,file_name=mesh_file)
                    
                    
            label = label+1
                 
        for structure in sb_structures:
                
            # find the first instance of segmentation
            #segfiles = glob(delineation_dir+'/'+subject+'*_mask-'+structure+'_*.nii.gz')
            segfiles = glob(delineation_dir+'/'+subject+'_ses-1_acq-sb*_mask-'+structure+'*.nii.gz')
                
            seg = None
            usewb=False
            if len(segfiles)>0:
                seg = load_volume(segfiles[0])
            else:
                segfiles = glob(delineation_dir+'/'+subject+'_ses-1_acq-wb*_mask-'+structure+'*.nii.gz')
                if len(segfiles)>0:
                    seg = load_volume(segfiles[0])
                    usewb=True
                
            if seg is not None:    
                print('load :'+segfiles[0])
                
                seg = nibabel.as_closest_canonical(seg)
                #nighres.io.save_volume(labeling_dir+'/'+subject+'_reorient_mask-'+structure+'.nii.gz',seg)

                if usewb:
                    seg = nighres.io.load_volume(nighres.registration.apply_coordinate_mappings(image=seg,
                                    mapping1=map_wb, interpolation='linear',
                                    save_data=True, file_name=subject+'_mask-'+structure+'_def-qmri2_slab-atlas.nii.gz',
                                    output_dir=atlas_dir)['result'])
                else:
                    seg = nighres.io.load_volume(nighres.registration.apply_coordinate_mappings(image=seg,
                                    mapping1=map_sb, interpolation='linear',
                                    save_data=True, file_name=subject+'_mask-'+structure+'_def-qmri2_slab-atlas.nii.gz',
                                    output_dir=atlas_dir)['result'])
                
                lvlset = nighres.io.load_volume(nighres.surface.probability_to_levelset(probability_image=seg,
                                        save_data=True, file_name=subject+'_mask-'+structure+'_def-qmri2_slab-atlas.nii.gz',
                                        output_dir=atlas_dir)['result']).get_fdata()
                
                if label_map is None:
                    min_dist = lvlset
                    label_map = nibabel.Nifti1Image((lvlset<=0)*label, seg.affine, seg.header)
                    min_dist[(lvlset>0)] = 0
                else:
                    data = numpy.maximum(label_map.get_fdata(),(lvlset<min_dist)*label)
                    label_map = nibabel.Nifti1Image(data, seg.affine, seg.header)
                    min_dist = numpy.minimum(min_dist,lvlset)
                
                if build_mesh:
                    mesh_file = subject+'_mesh-'+structure+'.vtk'
                    mesh = nighres.surface.levelset_to_mesh(lvlset,save_data=True,output_dir=mesh_dir,file_name=mesh_file)
                    
                    
            label = label+1
                 
        save_volume(label_file, label_map)
        label_map = None