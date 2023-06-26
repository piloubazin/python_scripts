import nighres
import numpy
import nibabel
from glob import glob
import os
import subprocess
import math
import scipy.ndimage
import ants

main_dir='/home/Public/jpnd/'
    
#in_dir=main_dir+'data/leipzig/traveling_phantom_7T/mpm_result/'
in_dir=main_dir+'data/leipzig/15636.fa/20230331/mpm/0p7/Results/'
ext_dir=main_dir+'data/leipzig/15636.fa/20230331/mpm/0p7/Results/Supplementary/'
ref_dir=main_dir+'data/amsterdam/traveling_phantom_7T/'
proc_dir=main_dir+'processing/15636.fa/mpm0p7/'

atlas_dir = '/home/pilou/Projects/Massp-Priors/'


# MPM input names
mpm_qr1_file = in_dir+'meas_MID00032_pdw_kp_mtflash3d_v1s_0p7_20230331-111133_rec-loraksRsos_echo-1_mt-off_part-mag_R1.nii.gz'
mpm_qr2s_file = in_dir+'meas_MID00032_pdw_kp_mtflash3d_v1s_0p7_20230331-111133_rec-loraksRsos_echo-1_mt-off_part-mag_R2s_WLS1.nii.gz'
mpm_qpd_file = in_dir+'meas_MID00032_pdw_kp_mtflash3d_v1s_0p7_20230331-111133_rec-loraksRsos_echo-1_mt-off_part-mag_PD.nii.gz'

# helper files?
mpm_T1w_file = ext_dir+'meas_MID00032_pdw_kp_mtflash3d_v1s_0p7_20230331-111133_rec-loraksRsos_echo-1_mt-off_part-mag_T1w_WLS1fit_TEzero.nii.gz'
mpm_PDw_file = ext_dir+'meas_MID00032_pdw_kp_mtflash3d_v1s_0p7_20230331-111133_rec-loraksRsos_echo-1_mt-off_part-mag_PDw_WLS1fit_TEzero.nii.gz'


# main processing steps

if (not os.path.exists(proc_dir)):
    os.makedirs(proc_dir)

print("\n0. MPM clean-up\n")

# N4 correction first on T1w, PDw images

# build a custom mask for N4, it seems the standard one crops data
n4mask_file = mpm_PDw_file.replace(ext_dir,proc_dir).replace('.nii.gz','_n4mask.nii.gz')
if (not os.path.isfile(n4mask_file)):
    n4mask = nighres.io.load_volume(mpm_PDw_file)
    n4mask = nibabel.Nifti1Image(n4mask.get_fdata()>0.01*numpy.max(n4mask.get_fdata()),n4mask.affine, n4mask.header)
    nighres.io.save_volume(n4mask_file, n4mask)

t1n4_file = mpm_T1w_file.replace(ext_dir,proc_dir).replace('.nii.gz','_n4.nii.gz')
if (not os.path.isfile(t1n4_file)):
    img = ants.image_read(mpm_T1w_file)
    t1_n4 = ants.n4_bias_field_correction(img, mask=n4mask_file)
    ants.image_write(t1_n4, t1n4_file)

pdn4_file = mpm_PDw_file.replace(ext_dir,proc_dir).replace('.nii.gz','_n4.nii.gz')
if (not os.path.isfile(pdn4_file)):
    img = ants.image_read(mpm_PDw_file)
    pd_n4 = ants.n4_bias_field_correction(img, mask=n4mask_file)
    ants.image_write(pd_n4, pdn4_file)

# skullstrip the MPMs + remove negative values
    
skull = nighres.brain.intensity_based_skullstripping(main_image=t1n4_file, extra_image=pdn4_file,
                    noise_model='exponential', 
                    skip_zero_values=True,
                    iterate=False, dilate_mask=0, topology_lut_dir=None,
                    save_data=True, overwrite=False, output_dir=proc_dir)

brainmask = nighres.io.load_volume(skull['brain_mask']).get_fdata()

mpm_r1strip_file = proc_dir+'mpm0p7_r1_brain.nii.gz'
if (not os.path.isfile(mpm_r1strip_file)):
    print("Mask qR1")
    r1 = nighres.io.load_volume(mpm_qr1_file)
    r1strip = nibabel.Nifti1Image(numpy.maximum(0.0,numpy.minimum(2.0,brainmask*r1.get_fdata())), r1.affine, r1.header)
    r1strip = nibabel.as_closest_canonical(r1strip)
    nighres.io.save_volume(mpm_r1strip_file, r1strip)

mpm_r2strip_file = proc_dir+'mpm0p7_r2s_brain.nii.gz'
if (not os.path.isfile(mpm_r2strip_file)):
    print("Mask qR2s")
    r2 = nighres.io.load_volume(mpm_qr2s_file)
    r2strip = nibabel.Nifti1Image(numpy.maximum(0.0,numpy.minimum(200.0,brainmask*r2.get_fdata())), r2.affine, r2.header)
    r2strip = nibabel.as_closest_canonical(r2strip)
    nighres.io.save_volume(mpm_r2strip_file, r2strip)

mpm_pdstrip_file = proc_dir+'mpm0p7_pd_brain.nii.gz'
if (not os.path.isfile(mpm_pdstrip_file)):
    print("Mask qPD")
    pd = nighres.io.load_volume(mpm_qpd_file)
    pd_data = pd.get_fdata()
    pd_data[pd_data==0] = 200.0
    pd_data[pd_data<0] = 200.0
    pdstrip = nibabel.Nifti1Image(brainmask*pd_data, pd.affine, pd.header)
    pdstrip = nibabel.as_closest_canonical(pdstrip)
    nighres.io.save_volume(mpm_pdstrip_file, pdstrip)


print("\n2. MASSP for MPM\n")

template = nighres.data.download_AHEAD_template()

# two-step ANTs registration for maximum precision
ants1 = nighres.registration.embedded_antspy_multi(
                            source_images=[mpm_r1strip_file,mpm_r2strip_file],
                            target_images=[template['qr1'],template['qr2s']],
                            run_rigid=True, run_affine=True, run_syn=True,
                            rigid_iterations=10000,
                            affine_iterations=2000,
                            coarse_iterations=180, 
                            medium_iterations=60, fine_iterations=30,
                            cost_function='MutualInformation', 
                            interpolation='Linear',
                            regularization='High',
                            smooth_mask=0.33,
                            ignore_affine=True, 
                            save_data=True, file_name='smpm0p7_coreg2atlas1',
                            output_dir=proc_dir)
    
ants2 = nighres.registration.embedded_antspy_multi(
                            source_images=ants1['transformed_sources'],
                            target_images=[template['qr1'],template['qr2s']],
                            run_rigid=False, run_affine=False, run_syn=True,
                            coarse_iterations=180, 
                            medium_iterations=60, fine_iterations=30,
                            cost_function='MutualInformation', 
                            interpolation='Linear',
                            regularization='Medium',
                            mask_zero=False,
                            smooth_mask=0.33,
                            ignore_affine=True, 
                            save_data=True, file_name='mpm0p7_coreg2atlas2',
                            output_dir=proc_dir)

mapping = nighres.registration.apply_coordinate_mappings(ants1['mapping'],mapping1=ants2['mapping'],
                            save_data=True, overwrite=False, file_name='mpm0p6_coreg2atlas-mapping',
                            output_dir=proc_dir)

inverse = nighres.registration.apply_coordinate_mappings(ants2['inverse'],mapping1=ants1['inverse'],
                            save_data=True, overwrite=False, file_name='mpm0p6_coreg2atlas-inverse',
                            output_dir=proc_dir)

massp = nighres.parcellation.massp(target_images=[mpm_r1strip_file,mpm_r2strip_file,mpm_pdstrip_file],
                            map_to_target=inverse['result'], 
                            shape_atlas_probas=atlas_dir+'atlas10_31struct_qmri2fcm_massp-sproba.nii.gz', 
                            shape_atlas_labels=atlas_dir+'atlas10_31struct_qmri2fcm_massp-slabel.nii.gz',
                            skeleton_atlas_probas=atlas_dir+'atlas10_31struct_qmri2fcm_massp-kproba.nii.gz', 
                            skeleton_atlas_labels=atlas_dir+'atlas10_31struct_qmri2fcm_massp-klabel.nii.gz',
                            intensity_atlas_hist=atlas_dir+'atlas10_31struct_mpm_massp-chist_cspmax-chist.nii.gz',
                            max_iterations=120, max_difference=0.1,
                            intensity_prior=0.25,
                            save_data=True, file_name='mpm0p7',
                            output_dir=proc_dir, overwrite=False)
