import nighres
import numpy
import math
import nibabel
from PIL import Image
from nighres.io import load_volume, save_volume
from nighres.utils import _fname_4saving
import scipy.ndimage
import os
from nibabel import processing
import subprocess
import shutil

subject='Ahead_brain_122017'

# input data: 3D image in the blockface space
input_img = 'Ahead_brain_122017_blockface-image.nii.gz'

stain_num = 308
# important: the slice number should always be specified in terms of stain number
# (the translation to blockface number is done internally here)

stain_ref = 'Ahead_brain_122017_Silver_308.tif'
# the stain reference is needed to resize the new image to same dimensions as
# original stain in the resampling step

# some stains have been scanned in flipped L-R orientation, set to True to correct
correct_orientation = True

# output data: 2D image with stain dimensions and corresponding slice number
out_dir = 'blockface_to_stain/'

# resampling and cropping parameters (do not change)
stain_scale=7
scaling = (0.15,0.15)

# location of the transformation files
mapping_dir = 'mappings/'
mapping_base = 'Ahead_brain_152017_stain2bf-'
mapping_suffix = '_invmap.nii.gz'

# get the input image: here we assume a color TIFF file
if not os.path.isfile(input_img):
    print('file not found')
else:
    mapping = mapping_dir+mapping_base+str(stain_num).zfill(3)+mapping_suffix
    if not os.path.isfile(mapping):
        print('mapping '+mapping+' not found')
    else:
        # corresponding blockface slice
        bfnum = stain_num+3
        if (stain_num<124): bfnum = stain_num + 16
        elif (stain_num<440): bfnum = stain_num + 17
        else: bfnum = stain_num + 16
        
        # extract 2D slice
        bf_img = load_volume(input_img)
        bf_slice = nibabel.Nifti1Image(bf_img.get_fdata()[:,:,bfnum],
                                       bf_img.affine, bf_img.header)
        
        inv = nighres.registration.apply_coordinate_mappings_2d(image=bf_slice, 
                                                          mapping1=mapping,
                                                          interpolation="nearest", 
                                                          padding="zero",
                                                          save_data=False)
    
        image = inv['result'].get_fdata()
            
        # flip if needed
        if not correct_orientation: image = numpy.flip(image, axis=0)    
        
        # scale up by stain_scale
        slice_expand = numpy.zeros((image.shape[0]*stain_scale,image.shape[1]*stain_scale))
        for dx in range(stain_scale):
            for dy in range(stain_scale):
                slice_expand[dx:slice_expand.shape[0]:stain_scale,
                             dy:slice_expand.shape[1]:stain_scale] = image

        # open original stain to get dimensions
        # note: TIFF files have X and Y dims transposed wrt Nifti files
        ref_img = numpy.array(Image.open(stain_ref).convert(mode='RGB'))
        print('resize: '+str(image.shape)+' -> '+str(slice_expand.shape)+' / '+str(ref_img.shape))
        slice_expand = slice_expand[0:ref_img.shape[1],0:ref_img.shape[0]]
            
        # save the result as nifti image
        header = nibabel.Nifti1Header()
        header.set_data_shape(slice_expand.shape)
        header.set_zooms(scaling)
            
        affine = numpy.eye(4)
        affine[0,3] = -slice_expand.shape[0]/2.0
        affine[1,3] = -slice_expand.shape[1]/2.0
            
        output_nifti = nibabel.Nifti1Image(slice_expand,affine=affine,header=header)
        out_file = out_dir+_fname_4saving(rootfile=input_img, 
                                    suffix='map2stain-'+str(stain_num).zfill(3), 
                                    ext='nii.gz')
        print('save to: '+out_file)
        
        save_volume(out_file, output_nifti)
