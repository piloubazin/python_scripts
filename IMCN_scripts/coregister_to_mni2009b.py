
import nighres
import numpy
import math
import nibabel
import os
import shutil

subject='Ahead_brain_122017'

blockface = 'Ahead_brain_122017_blockface-image.nii.gz'

# here we assume to have a copy of the MNI 2009B 0.5mm template
# as well as a scaled brain mask derived from 2009A 1mm
mni2009b_mask = 'mni_icbm152_t1_tal_nlin_asym_09a_mask_scaled09b.nii.gz'
mni2009b_t1 = 'mni_icbm152_nlin_asym_09b/mni_icbm152_t1_tal_nlin_asym_09b_hires.nii.gz'

silver = 'Ahead_brain_122017_Bielschowsky-interpolated.nii.gz'
thio = 'Ahead_brain_122017_thionin-interpolated.nii.gz'
parv = 'Ahead_brain_122017_parvalbumin-interpolated.nii.gz'
calb = 'Ahead_brain_122017_calbindin-interpolated.nii.gz'
calret = 'Ahead_brain_122017_calretinin-interpolated.nii.gz'

qr1 = 'Ahead_brain_122017_MRI-quantitative-R1.nii.gz'
qr2s = 'Ahead_brain_122017_MRI-quantitative-R2star.nii.gz'
qpd = 'Ahead_brain_122017_MRI-proton-density.nii.gz'


output_dir = 'mni2009b/'
  
# mask the MNI template  
mni2009b_t1masked = output_dir+'mni_icbm152_t1_tal_nlin_asym_09b_hires_masked.nii.gz'
if (os.path.isfile(mni2009b_t1masked)):
        print('MNI template masking: done')
else:
    t1 = nibabel.funcs.as_closest_canonical(nighres.io.load_volume(mni2009b_t1))
    mask = nibabel.funcs.as_closest_canonical(nighres.io.load_volume(mni2009b_mask)).get_fdata()
    
    t1_masked = nibabel.Nifti1Image(mask*t1.get_fdata(),affine=t1.affine,header=t1.header)
    nighres.io.save_volume(mni2009b_t1masked,t1_masked)

# co-register blockface to MNI template
output = output_dir+subject+'_bf2mni2009b.nii.gz'

mni = nighres.registration.embedded_antsreg(source_image=blockface, 
                    target_image=mni2009b_t1masked,
                    run_rigid=True,
                    rigid_iterations=1000,
                    run_affine=True,
                    affine_iterations=2000,
                    run_syn=True,
                    coarse_iterations=200,
                    medium_iterations=100, fine_iterations=50,
					cost_function='MutualInformation',
					interpolation='NearestNeighbor',
					regularization='High',
					convergence=1e-6,
					mask_zero=False,
					ignore_affine=True,
                    save_data=True, overwrite=False,
                    file_name=output)

# Transform histology images
nighres.registration.apply_coordinate_mappings(silver, mni['mapping'], interpolation="nearest", padding="closest",
                        save_data=True, overwrite=False, output_dir=output_dir, file_name=subject+'_silver')

nighres.registration.apply_coordinate_mappings(thio, mni['mapping'], interpolation="nearest", padding="closest",
                        save_data=True, overwrite=False, output_dir=output_dir, file_name=subject+'_thio')

nighres.registration.apply_coordinate_mappings(parv, mni['mapping'], interpolation="nearest", padding="closest",
                        save_data=True, overwrite=False, output_dir=output_dir, file_name=subject+'_parv')

nighres.registration.apply_coordinate_mappings(calb, mni['mapping'], interpolation="nearest", padding="closest",
                        save_data=True, overwrite=False, output_dir=output_dir, file_name=subject+'_calb')

nighres.registration.apply_coordinate_mappings(calret, mni['mapping'], interpolation="nearest", padding="closest",
                        save_data=True, overwrite=False, output_dir=output_dir, file_name=subject+'_calret')

nighres.registration.apply_coordinate_mappings(qr1, mni['mapping'], interpolation="nearest", padding="closest",
                        save_data=True, overwrite=True, output_dir=output_dir, file_name=subject+'_qr1')

nighres.registration.apply_coordinate_mappings(qr2s, mni['mapping'], interpolation="nearest", padding="closest",
                        save_data=True, overwrite=True, output_dir=output_dir, file_name=subject+'_qr2s')

nighres.registration.apply_coordinate_mappings(qpd, mni['mapping'], interpolation="nearest", padding="closest",
                        save_data=True, overwrite=True, output_dir=output_dir, file_name=subject+'_qpd')
