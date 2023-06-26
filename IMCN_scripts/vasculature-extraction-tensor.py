import nighres
import numpy
import nibabel
from glob import glob
import os
import subprocess
import scipy.ndimage
mrtrix = '/home/pilou/Software/mrtrix3/bin/'

main_dir='/home/pilou/Projects/Ahead-Database/Vascular-Modeling/'
in_dir='/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/'
out_dir=main_dir+'process_qmri2/'
tens_dir=main_dir+'process_tensor/'
final_dir=main_dir+'vessel-maps/'


subjects = range(91,92)

for num in subjects:
    subject = 'sub-'+str(num).zfill(3)

    print("process subject "+subject)
    
    r1_file = in_dir+subject+'/ses-1/anat/wb/qmri/'+subject+'_ses-1_acq-wb_mod-r1hz_orient-std_brain.nii.gz'
    r2_file = in_dir+subject+'/ses-1/anat/wb/qmri/'+subject+'_ses-1_acq-wb_mod-r2hz_orient-std_brain.nii.gz'
    qsm_file = in_dir+subject+'/ses-1/anat/wb/qmri/'+subject+'_ses-1_acq-wb_mod-qsm_orient-std_brain.nii.gz'

    if (os.path.isfile(r1_file) and os.path.isfile(r2_file) and os.path.isfile(qsm_file)):

        print("Mask qR1")
        r1 = nighres.io.load_volume(r1_file)
        mask = scipy.ndimage.binary_erosion((r1.get_fdata()>0), iterations=5)
        mask_img = nibabel.Nifti1Image(mask, r1.affine, r1.header)
        
        print("1. Vasculature from R2*")
        
        vasc1 = nighres.filtering.multiscale_vessel_filter(r2_file, structure_intensity='bright', filterType = 'RRF',
                                        propagationtype = 'diffusion', threshold = 0.5,
                                        factor = 0.5,max_diff = 0.001,max_itr = 100,
                                        scale_step = 1.0,
                                        scales = 3,
                                        prior_image=mask_img,
                                        invert_prior=False,
                                        save_data=True,
                                        overwrite=False,
                                        output_dir=out_dir)   
        
        print("2. Vasculature from R1")
        
        vasc2 = nighres.filtering.multiscale_vessel_filter(r1_file, structure_intensity='bright', filterType = 'RRF',
                                        propagationtype = 'diffusion', threshold = 0.5,
                                        factor = 0.5,max_diff = 0.001,max_itr = 100,
                                        scale_step = 1.0,
                                        scales = 3,
                                        prior_image=mask_img,
                                        invert_prior=False,
                                        save_data=True,
                                        overwrite=False,
                                        output_dir=out_dir)   
        
        print("3. QSM")
            
        vasc3 = nighres.filtering.multiscale_vessel_filter(qsm_file, structure_intensity='bright', filterType = 'RRF',
                                        propagationtype = 'diffusion', threshold = 0.5,
                                        factor = 0.5,max_diff = 0.001,max_itr = 100,
                                        scale_step = 1.0,
                                        scales = 3,
                                        prior_image=mask_img,
                                        invert_prior=False,
                                        save_data=True,
                                        overwrite=False,
                                        output_dir=out_dir)   

        print("4. Structure tensors")

        overwrite=True

        tensor_file = r2_file.replace(in_dir+subject+'/ses-1/anat/wb/qmri/',tens_dir).replace('.nii.gz','_st.nii.gz')
        if ((os.path.exists(tensor_file)) & (overwrite == False)):
            print('\n Generate tensor: done (set overwrite to True to recompute)')
        else:
            print('\n Generate tensor')
            os.makedirs(tens_dir, exist_ok=True)
            img = nighres.io.load_volume(r2_file)
            data = img.get_fdata()
            
            dx = 0.5*(numpy.roll(data,1,axis=0) - numpy.roll(data,-1,axis=0))
            dy = 0.5*(numpy.roll(data,1,axis=1) - numpy.roll(data,-1,axis=1))
            dz = 0.5*(numpy.roll(data,1,axis=2) - numpy.roll(data,-1,axis=2))
            
            dxx = numpy.roll(data,1,axis=0) - 2.0*data + numpy.roll(data,-1,axis=0)
            dyy = numpy.roll(data,1,axis=1) - 2.0*data + numpy.roll(data,-1,axis=1)
            dzz = numpy.roll(data,1,axis=2) - 2.0*data + numpy.roll(data,-1,axis=2)
            
            dxy = 0.25*(numpy.roll(data,[1,1],axis=[0,1]) - numpy.roll(data,[-1,1],axis=[0,1]) - numpy.roll(data,[1,-1],axis=[0,1]) + numpy.roll(data,[-1,-1],axis=[0,1]))
            dxz = 0.25*(numpy.roll(data,[1,1],axis=[0,2]) - numpy.roll(data,[-1,1],axis=[0,2]) - numpy.roll(data,[1,-1],axis=[0,2]) + numpy.roll(data,[-1,-1],axis=[0,2]))
            dyz = 0.25*(numpy.roll(data,[1,1],axis=[1,2]) - numpy.roll(data,[-1,1],axis=[1,2]) - numpy.roll(data,[1,-1],axis=[1,2]) + numpy.roll(data,[-1,-1],axis=[1,2]))
            
            # tensor is built from shape operator (dx,dy,dz)^-1(dxx,dyy,dzz), requires some closed form inversion
            norm = numpy.sqrt(dx*dx+dy*dy+dz*dz)
            
            tensor = numpy.zeros((data.shape[0],data.shape[1],data.shape[2],6))
            tensor[:,:,:,0] = norm-dx*dx
            tensor[:,:,:,1] = norm-dy*dy
            tensor[:,:,:,2] = norm-dz*dz
            tensor[:,:,:,3] = dxy*dxy
            tensor[:,:,:,4] = dxz*dxz
            tensor[:,:,:,5] = dyz*dyz
            
            tensor_img = nibabel.Nifti1Image(tensor, img.affine, img.header)
            nighres.io.save_volume(tensor_file, tensor_img)
        
        overwrite=True
        
        fa_file = tensor_file.replace('.nii.gz','-fa.nii.gz')
        adc_file = tensor_file.replace('.nii.gz','-adc.nii.gz')
        vec_file = tensor_file.replace('.nii.gz','-vec.nii.gz')
        if ((os.path.exists(fa_file)) & (os.path.exists(adc_file)) & (os.path.exists(vec_file)) & (overwrite == False)):
            print('\n Generate metrics: done (set overwrite to True to recompute)')
        else:
            print('\n Generate metrics')
     
            command = mrtrix+'tensor2metric -force -fa '+fa_file+' -adc '+adc_file+' -vec '+vec_file+' '+tensor_file
                    
            print(command+'\n')
            try:
                subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
                raise subprocess.CalledProcessError(msg)
                
        overwrite=False

        tck_file = tensor_file.replace('.nii.gz','.tck')
        if ((os.path.exists(tck_file)) & (overwrite == False)):
            print('\n Generate tracts: done (set overwrite to True to recompute)')
        else:
            print('\n Generate tracts')
     
            command = mrtrix+'tckgen '+tensor_file+' '+tck_file+' -seed_image '+vasc1['segmentation']\
                            +' -nthreads 4 -algorithm FACT -select 10000 -force -cutoff 0.05'
             
            print(command+'\n')
            try:
                subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
                raise subprocess.CalledProcessError(msg)
                
        map_file = tensor_file.replace('.nii.gz','-map.nii.gz')
        if ((os.path.exists(map_file)) & (overwrite == False)):
            print('\n Generate maps: done (set overwrite to True to recompute)')
        else:
            print('\n Generate maps')
     
            command = mrtrix+'tckmap -contrast scalar_map -image '+vasc1['probability']+' -force '\
                          +'-template '+vasc1['probability']+' '+tck_file+' '+map_file

            print(command+'\n')
            try:
                subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
                raise subprocess.CalledProcessError(msg)
                