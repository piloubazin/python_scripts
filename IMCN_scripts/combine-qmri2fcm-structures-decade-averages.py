import nibabel
import nighres
import glob
import numpy
import os
from nighres.io import load_volume, save_volume
import sys

#structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
#              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r']
structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r',
              'gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r','tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r',
              'vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r','ic_hem-l','ic_hem-r','vta_hem-l','vta_hem-r',
              'fx_hem-lr','pag_hem-l','pag_hem-r','ppn_hem-l','ppn_hem-r','cl_hem-l','cl_hem-r']


ahead_brain = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/ahead_qmri2_mni09b_med_r1map_n105.nii.gz'

declb = ['18to30','31to40','41to50','51to60','61to70','71to80']
declb = ['18to30','31to40']
declb = ['41to50','51to60','61to70','71to80']

declb = ['61to80','41to60','18to40','18to80']

header = load_volume(ahead_brain).header
affine = load_volume(ahead_brain).affine
ahead_geom = (header.get_data_shape()[0],header.get_data_shape()[1],header.get_data_shape()[2])

for decade in declb:
    maxproba_name = 'ahead-qmri2fcm_avg-maxproba_decade-'+decade+'.nii.gz'
    maxlabel_name = 'ahead-qmri2fcm_avg-maxlabel_decade-'+decade+'.nii.gz'
    seglabel_name = 'ahead-qmri2fcm_avg-bestlabel_decade-'+decade+'.nii.gz'

    if not (os.path.exists(maxproba_name) and os.path.exists(maxlabel_name) and os.path.exists(seglabel_name) ):
        print("combining all atlases for decade "+decade)
        
        maxproba = numpy.zeros(ahead_geom)
        maxlabel = numpy.zeros(ahead_geom)
        seglabel = numpy.zeros(ahead_geom)
        
        background = glob.glob('/home/pilou/Projects/Ahead-Database/Automated-Parcellation/ahead-qmri2fcm_avg-background_decade-'+decade+'*.nii.gz')[0]    
        background = load_volume(background).get_data()
        
        print("averaging")
        for index,structure in enumerate(structures):
            print("structure: "+structure)
            
            atlas = glob.glob('/home/pilou/Projects/Ahead-Database/Automated-Parcellation/ahead-qmri2fcm_avg-'+structure+'_decade-'+decade+'*.nii.gz')[0]    
                    
            img = load_volume(atlas)
            data = img.get_data()
        
            seglabel = (data>background)*( (data>maxproba)*(index+1) + (data<=maxproba)*seglabel) + (data<=background)*seglabel
            maxlabel = (data>maxproba)*(index+1) + (data<=maxproba)*maxlabel
            maxproba = (data>maxproba)*data + (data<=maxproba)*maxproba
            
        maxproba_img = nibabel.Nifti1Image(maxproba, affine, header)
        save_volume(maxproba_name, maxproba_img)
            
        maxlabel_img = nibabel.Nifti1Image(maxlabel, affine, header)
        save_volume(maxlabel_name, maxlabel_img)
        
        seglabel_img = nibabel.Nifti1Image(seglabel, affine, header)
        save_volume(seglabel_name, seglabel_img)
    
    maxproba_name = 'ahead-qmri2fcm_std-maxproba_decade-'+decade+'.nii.gz'
    if not os.path.exists(maxproba_name):
        print("standard deviation")
        maxproba = numpy.zeros(ahead_geom)
        for index,structure in enumerate(structures):
            print("structure: "+structure)
            
            atlas = glob.glob('/home/pilou/Projects/Ahead-Database/Automated-Parcellation/ahead-qmri2fcm_std-'+structure+'_decade-'+decade+'*.nii.gz')[0]    
                    
            img = load_volume(atlas)
            data = img.get_data()
        
            #maxproba = (data>maxproba)*data + (data<=maxproba)*maxproba
            maxproba = (maxlabel==index+1)*data + (maxlabel!=index+1)*maxproba
            
        maxproba_img = nibabel.Nifti1Image(maxproba, affine, header)
        save_volume(maxproba_name, maxproba_img)


    nighres.surface.parcellation_to_meshes(seglabel_name, connectivity="18/6", 
                     spacing = 0.0, smoothing=1.0,
                     save_data=True, overwrite=False,
                     file_name='ahead-qmri2fcm_avg-bestlabel_decade-'+decade)

