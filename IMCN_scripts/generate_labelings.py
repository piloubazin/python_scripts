import nibabel
import nighres
import glob
import numpy
import os
from nighres.io import load_volume, save_volume

delineation_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/manual_delineations/best'
labeling_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/label_maps'
mesh_dir = '/home/pilou/Projects/Subcortex/MultiAtlas-Segmentation/meshes'

#subjects = ['sub-000','sub-001','sub-002','sub-004','sub-008','sub-012','sub-024','sub-031','sub-033','sub-040']
#structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r']

#subjects = ['sub-000','sub-002','sub-004','sub-008','sub-012']
#structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
#              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4',
#              'vta_hem-l','vta_hem-r','scg_hem-l','scg_hem-r','amg_hem-l','amg_hem-r','pag_hem-l','pag_hem-r','lh_hem-l','lh_hem-r','fx']

subjects = ['sub-000','sub-001','sub-002','sub-004','sub-008','sub-012','sub-024','sub-031','sub-033','sub-040']
#subjects = ['sub-033']

#structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
#              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r',
#              'ico_hem-l','ico_hem-r','cl_hem-l','cl_hem-r','fx_hem-lr','lh_hem-l','lh_hem-r','pag_hem-l','pag_hem-r','ppn_hem-l','ppn_hem-r',
#              'scg_hem-l','scg_hem-r','sco_hem-l','sco_hem-r','vta_hem-l','vta_hem-r']
#structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
#              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4']

#structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
#              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r',
#              'ico_hem-l','ico_hem-r','cl_hem-l','cl_hem-r','fx_hem-lr','pag_hem-l','pag_hem-r','sco_hem-l','sco_hem-r','vta_hem-l','vta_hem-r']

#structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
#              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r',
#              'ic_hem-l','ic_hem-r','vta_hem-l','vta_hem-r','fx_hem-lr','ppn_hem-l','ppn_hem-r','cl_hem-l','cl_hem-r','lh_hem-l','lh_hem-r',
#              'ico_hem-l','ico_hem-r','pag_hem-l','pag_hem-r','sco_hem-l','sco_hem-r','scg_hem-l','scg_hem-r']

#structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
#              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r',
#              'ic_hem-l','ic_hem-r','vta_hem-l','vta_hem-r','fx_hem-lr','ppn_hem-l','ppn_hem-r','pag_hem-l','pag_hem-r',
#              'ico_hem-l','ico_hem-r','sco_hem-l','sco_hem-r','cl_hem-l','cl_hem-r']

#structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
#              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r',
#              'ic_hem-l','ic_hem-r','vta_hem-l','vta_hem-r','fx_hem-lr','pag_hem-l','pag_hem-r','ppn_hem-l','ppn_hem-r',
#              'cl_hem-l','cl_hem-r']

#structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
#              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r',
#              'ic_hem-l','ic_hem-r','vta_hem-l','vta_hem-r','fx_hem-lr','pag_hem-l','pag_hem-r','ppn_hem-l','ppn_hem-r',
#              'cl_hem-l','cl_hem-r','cow_hem-lr']

structures = ['str_hem-l','str_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r',
              'ic_hem-l','ic_hem-r',
              'cl_hem-l','cl_hem-r','cow_hem-lr']

build_mesh = False

label_map = None
min_dist = None
for subject in subjects:
    
    label_file = labeling_dir+'/'+subject+'_labeling-'+str(len(structures))+'.nii.gz'
    
    label = 1
    for structure in structures:
            
        # find the first instance of segmentation
        segfiles = glob.glob(delineation_dir+'/'+subject+'*_mask-'+structure+'_*.nii.gz')
        
        seg = None
        if len(segfiles)>0:
            seg = load_volume(segfiles[0])
            
        if seg is not None:    
            print('load :'+segfiles[0])
            data = seg.get_data()
            
            lvlset = nighres.surface.probability_to_levelset(seg)['result'].get_data()
            
            if label_map is None:
                min_dist = lvlset
                label_map = nibabel.Nifti1Image((lvlset<=0)*label, seg.affine, seg.header)
                min_dist[(lvlset>0)] = 0
            else:
                data = numpy.maximum(label_map.get_data(),(lvlset<min_dist)*label)
                label_map = nibabel.Nifti1Image(data, seg.affine, seg.header)
                min_dist = numpy.minimum(min_dist,lvlset)
            
            if build_mesh:
                mesh_file = subject+'_mesh-'+structure+'.vtk'
                mesh = nighres.surface.levelset_to_mesh(lvlset,save_data=True,output_dir=mesh_dir,file_name=mesh_file)
                
                
        label = label+1
             
    save_volume(label_file, label_map)
    label_map = None