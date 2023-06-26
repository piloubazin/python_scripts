import nibabel as nb
import nighres as nr
import glob
import numpy
import os
from nighres.io import load_volume, save_volume


# List of subject IDs to be processed, assuming they all come from the Subcortex database (BIDS version)
#subjects = sorted(glob.glob('./r1maps/*.nii.gz'))
#mni_brain = '/usr/share/fsl/data/standard/MNI152_T1_0.5mm.nii.gz'            
#mni_brain = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/MNI152_T1_0.5mm_brain_d4mm.nii.gz'            
mni_brain = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/mni_icbm152_t1_tal_nlin_asym_09b_hires_brain_d4mm.nii.gz'

tmp_dir = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/tmp/'

subsample = 4
maxsubject = 110

header = load_volume(mni_brain).header
affine = load_volume(mni_brain).affine
mni_geom = (header.get_data_shape()[0],header.get_data_shape()[1],header.get_data_shape()[2])

median = numpy.zeros(mni_geom)
iqr = numpy.zeros(mni_geom)

# first check for number of subjects
nsubjects = 0
for subject in range(0,maxsubject): 
    subject = str(subject).zfill(3)
    whole_brain = '/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/sub-'+subject\
    			+'/ses-1/anat/wb/qmri/sub-'+subject+'_ses-1_acq-wb_mod-r1hz_orient-std_brain.nii.gz'
    
    mapping = os.path.join('/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/mni2009b-qmri2-mappings/',
                            'sub-'+subject+'_ses-1_acq-wb_mod-r1hz_orient-std_brain_map-mni09b_ants-map.nii.gz')
    
    if os.path.exists(whole_brain) is False: print("error: whole brain doesn't exist")
    elif os.path.exists(mapping) is False: print("error: mapping doesn't exist")
    else:
        nsubjects = nsubjects+1


# create a subsampled data set because memory cannot handle the whole sized image X 105
for dx in range(subsample) :
    for dy in range(subsample) :
        for dz in range(subsample) :
            # coordinates to sample from
            xmin = int(dx*numpy.ceil(mni_geom[0]/subsample))
            xmax = int(numpy.minimum(mni_geom[0],(dx+1)*numpy.ceil(mni_geom[0]/subsample)))
            
            ymin = int(dy*numpy.ceil(mni_geom[1]/subsample))
            ymax = int(numpy.minimum(mni_geom[1],(dy+1)*numpy.ceil(mni_geom[1]/subsample)))
            
            zmin = int(dz*numpy.ceil(mni_geom[2]/subsample))
            zmax = int(numpy.minimum(mni_geom[2],(dz+1)*numpy.ceil(mni_geom[2]/subsample)))
            
            print('sub-image: ['+str(xmin)+','+str(xmax)+']x['+str(ymin)+','+str(ymax)+']x['+str(zmin)+','+str(zmax)+']')
            
            sub_geom = (xmax-xmin,ymax-ymin,zmax-zmin,nsubjects)
            histos = numpy.zeros(sub_geom)
            
            subnum=0
            for subject in range(0,maxsubject) : 
                subject = str(subject).zfill(3)
                whole_brain = '/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/sub-'+subject\
                            +'/ses-1/anat/wb/qmri/sub-'+subject+'_ses-1_acq-wb_mod-r1hz_orient-std_brain.nii.gz'
                
                mapping = os.path.join('/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/mni2009b-qmri2-mappings/',
                                        'sub-'+subject+'_ses-1_acq-wb_mod-r1hz_orient-std_brain_map-mni09b_ants-map.nii.gz')
    
                print("adding: "+whole_brain+" to atlas")
                if os.path.exists(whole_brain) is False: print("error: whole brain doesn't exist")
                elif os.path.exists(mapping) is False: print("error: mapping doesn't exist")
                else:
                    mapped = nr.registration.apply_coordinate_mappings(whole_brain, mapping, 
                                                    interpolation="nearest", padding="closest",
                                                    save_data=True,output_dir=tmp_dir)

                    histos[:,:,:,subnum] = nr.io.load_volume(mapped['result']).get_fdata()[xmin:xmax,ymin:ymax,zmin:zmax]
                    subnum = subnum+1
                
            # compute the interesting statistics
            median[xmin:xmax,ymin:ymax,zmin:zmax] = numpy.percentile(histos, 50, axis=3)
            iqr[xmin:xmax,ymin:ymax,zmin:zmax] = numpy.percentile(histos, 75, axis=3) \
                                                -numpy.percentile(histos, 25, axis=3)


med_img = nb.Nifti1Image(median, affine, header)
med_name = 'ahead_qmri2_mni09b_med_r1map_n'+str(nsubjects)+'.nii.gz'
save_volume(med_name, med_img)

iqr_img = nb.Nifti1Image(iqr, affine, header)
iqr_name = 'ahead_qmri2_mni09b_iqr_r1map_n'+str(nsubjects)+'.nii.gz'
save_volume(iqr_name, iqr_img)

