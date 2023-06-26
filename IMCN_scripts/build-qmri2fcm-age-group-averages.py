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
              
labels = [1,2,3,4,5,6,7,8,
          9,10,11,12,13,14,15,16,
          17,18,19,20,21,22,23,24,
          25,26,27,28,29,30,31]
            
decade1830 = [0,1,2,3,4,5,6,7,8,12,18,19,20,21,24,25,26,27,28,30,33,34,35,36,40,41,42,43,45,46,50,53,55,59,61,73,104,105,106,107,108,109]
decade3140 = [10,11,31,75,80,83,84,85,90,92,93,95]
decade4150 = [22,32,38,39,44,56,66,82,89,91,94,96,97]
decade5160 = [37,48,52,54,58,62,63,65,69,71,79,86]
decade6170 = [9,14,16,17,23,47,49,60,64,67,68,72,100]
decade7180 = [51,74,76,77,78,81,87,88,98,99,101,102,103]

decades = [decade1830,decade3140,decade4150,decade5160,decade6170,decade7180]
declb = ['18to30','31to40','41to50','51to60','61to70','71to80']

# doable, but overkill given the existing data, unless stdev is interesting...
decade1840 = decade1830+decade3140
decade4160 = decade4150+decade5160
decade6180 = decade6170+decade7180
decade1880 = decade1830+decade3140+decade4150+decade5160+decade6170+decade7180

decades = [decade1840,decade4160,decade6180,decade1880]
declb = ['18to40','41to60','61to80','18to80']

decades = [decade1880]
declb = ['18to80']

ahead_brain = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/ahead_qmri2_mni09b_med_r1map_n105.nii.gz'

header = load_volume(ahead_brain).header
affine = load_volume(ahead_brain).affine
ahead_geom = (header.get_data_shape()[0],header.get_data_shape()[1],header.get_data_shape()[2])

if (len(sys.argv)>1):
    print("atlasing "+sys.argv[1])

for dec,decade in enumerate(decades):
    print("decade "+declb[dec])
    for index,structure in enumerate(structures):
#        if (len(sys.argv)>1 and structure==sys.argv[1]) or (len(sys.argv)==1):
        if index>26:
            average = numpy.zeros(ahead_geom)
            nsubjects = 0
            for subject in decade: 
                subject = str(subject).zfill(3)
                parcellation = os.path.join('/home/pilou/Projects/Ahead-Database/Automated-Parcellation/qmri2fcm_all/',
                                        'sub-'+subject+'_ses-1_tim_massp-label.nii.gz')
            
                mapping = os.path.join('/home/pilou/Projects/Ahead-Database/Automated-Parcellation/qmri2fcm_proc/',
                                        'sub-'+subject+'_ses-1_reg2ahead-mapping_def-img.nii.gz')
            
                print("adding: "+parcellation+" to "+structure+" atlas")
                if os.path.exists(parcellation) is False: print("error: whole brain doesn't exist")
                elif os.path.exists(mapping) is False: print("error: mapping doesn't exist")
                else:
                    # create a mask first
                    img = load_volume(parcellation)
                    data = img.get_data()
                    struct = nb.Nifti1Image(data==labels[index],img.affine,img.header)
                
                    mapped = nr.registration.apply_coordinate_mappings(struct, mapping, 
                                                    interpolation="linear", padding="zero",
                                                    save_data=False)
        
                    average = average + mapped['result'].get_data()
                    nsubjects = nsubjects+1
                
            if nsubjects>0: average = average / nsubjects
        
            avg_img = nb.Nifti1Image(average, affine, header)
            avg_name = 'ahead-qmri2fcm_avg-'+structure+'_decade-'+declb[dec]+'_n'+str(nsubjects)+'.nii.gz'
            save_volume(avg_name, avg_img)

            stdev = numpy.zeros(ahead_geom)
            nsubjects = 0
            for subject in decade: 
                subject = str(subject).zfill(3)
                parcellation = os.path.join('/home/pilou/Projects/Ahead-Database/Automated-Parcellation/qmri2fcm_all/',
                                        'sub-'+subject+'_ses-1_tim_massp-label.nii.gz')
            
                mapping = os.path.join('/home/pilou/Projects/Ahead-Database/Automated-Parcellation/qmri2fcm_proc/',
                                        'sub-'+subject+'_ses-1_reg2ahead-mapping_def-img.nii.gz')
            
                print("adding: "+parcellation+" to "+structure+" atlas")
                if os.path.exists(parcellation) is False: print("error: whole brain doesn't exist")
                elif os.path.exists(mapping) is False: print("error: mapping doesn't exist")
                else:
                    # create a mask first
                    img = load_volume(parcellation)
                    data = img.get_data()
                    struct = nb.Nifti1Image(data==labels[index],img.affine,img.header)
                
                    mapped = nr.registration.apply_coordinate_mappings(struct, mapping, 
                                                    interpolation="linear", padding="zero",
                                                    save_data=False)
        
                    stdev = stdev + (mapped['result'].get_data() - average)*(mapped['result'].get_data() - average)
                    nsubjects = nsubjects+1
                
            if nsubjects>0: stdev = numpy.sqrt(stdev / nsubjects)
            
            std_img = nb.Nifti1Image(stdev, affine, header)
            std_name = 'ahead-qmri2fcm_std-'+structure+'_decade-'+declb[dec]+'_n'+str(nsubjects)+'.nii.gz'
            save_volume(std_name, std_img)
        
