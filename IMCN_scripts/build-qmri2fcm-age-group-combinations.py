import nibabel as nb
import nighres as nr
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
decnb = [42,12,13,12,13,13]

grplb = '18to80'

header = load_volume(ahead_brain).header
affine = load_volume(ahead_brain).affine
ahead_geom = (header.get_data_shape()[0],header.get_data_shape()[1],header.get_data_shape()[2])

for structure in structures:
    print("combining all atlases for structure "+structure)
    
    avgproba = numpy.zeros(ahead_geom)
    grpnb = numpy.sum(decnb)
     
    print("averaging")
    for index,decade in enumerate(declb):
        print("decade: "+declb)
        
        atlas = glob.glob('/home/pilou/Projects/Ahead-Database/Automated-Parcellation/atlas_maps/qmri2fcm/ahead-qmri2fcm_avg-'+structure+'_decade-'+decade+'*.nii.gz')[0]    
        img = load_volume(atlas)
        data = img.get_data()
    
        avgproba = avgproba + decnb*data
        
    avgproba = avgproba/grpnb
        
    avgproba_img = nb.Nifti1Image(proba, affine, header)
    avgproba_name = 'ahead-qmri2fcm_avg-'+structure+'_decades-'+grplb+'_n'+str(grpnb)+'.nii.gz'
    save_volume(avgproba_name, avgproba_img)
        
