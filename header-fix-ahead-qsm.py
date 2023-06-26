import os
import shutil
import glob
import nibabel
import numpy

in_folder = '/home/public/HumanAtlas_BIDS/'
out_folder = '/home/public/HumanAtlas_BIDS/'
in_path = '/ses-1/anat/wb/qmri/'
out_path = '/ses-1/anat/wb/qmri/'

os.makedirs(out_folder,exist_ok=True)

subjects = glob.glob(in_folder+'sub-*/')
for subject in subjects:
    subject = subject.replace(in_folder,'').replace('/','')
    print(subject)
    
    header = glob.glob(in_folder+subject+in_path+'*mod-qsm*')
    if len(header)>0:
        header = header[0]
            
        data_files = glob.glob(in_folder+subject+in_path+'*chi_*')
        for data in data_files:
            print(data+" header <- "+header)
            
            img = nibabel.load(data)
            hdr = nibabel.load(header)
            
            img = nibabel.Nifti1Image(img.get_fdata(),hdr.affine,hdr.header)
            #data = data.replace('.nii.gz','_hfx.nii.gz')
            nibabel.save(img, data)
                
        data_files = glob.glob(in_folder+subject+in_path+'*lfs_*')
        for data in data_files:
            print(data+" header <- "+header)
            
            img = nibabel.load(data)
            hdr = nibabel.load(header)
            
            img = nibabel.Nifti1Image(img.get_fdata(),hdr.affine,hdr.header)
            #data = data.replace('.nii.gz','_hfx.nii.gz')
            nibabel.save(img, data)
        
