import nighres
import numpy
import nibabel
from glob import glob
import os
import subprocess
import math
import scipy.ndimage

main_dir='/home/Public/jpnd/'
    
#in_dir=main_dir+'data/leipzig/traveling_phantom_7T/mpm_result/'
in_dir=main_dir+'data/leipzig/traveling_phantom_7T/mpm_denoised/mpm_wb_c110/'
ext_dir=main_dir+'data/leipzig/traveling_phantom_7T/mpm_calc/'
ref_dir=main_dir+'data/amsterdam/traveling_phantom_7T/'
proc_dir=main_dir+'processing/traveling_phantom_7T/'
out_dir=main_dir+'derivatives/traveling_phantom_7T/'

# MPM input names
#mpm_qr1_file = in_dir+'s2022-12-15_19-21-200050-00001-00288-1_R1.nii.gz'
#mpm_qr2s_file = in_dir+'s2022-12-15_19-21-200050-00001-00288-1_R2s_WLS1.nii.gz'
#mpm_qpd_file = in_dir+'s2022-12-15_19-21-200050-00001-00288-1_PD.nii.gz'
mpm_qr1_file = in_dir+'s2022-12-15_19-21-200050-00001-00288-1_c110_lcpca-den_R1.nii.gz'
mpm_qr2s_file = in_dir+'s2022-12-15_19-21-200050-00001-00288-1_c110_lcpca-den_R2s_WLS1.nii.gz'
mpm_qpd_file = in_dir+'s2022-12-15_19-21-200050-00001-00288-1_c110_lcpca-den_PD.nii.gz'

# helper files?
mpm_T1w_file = ext_dir+'s2022-12-15_19-21-200050-00001-00288-1_T1w_WLS1fit_TEzero.nii.gz'
mpm_PDw_file = ext_dir+'s2022-12-15_19-21-200050-00001-00288-1_PDw.nii.gz'

# AHEAD input names
ahead_qr1_file = ref_dir+'sub-091_ses-1_acq-wb2_mod-r1hz_orient-std_brain.nii.gz'
ahead_qr2s_file = ref_dir+'sub-091_ses-1_acq-wb2_mod-r2hz_orient-std_brain.nii.gz'
ahead_qsm_file = ref_dir+'sub-091_ses-1_acq-wb2_mod-qsm_orient-std_brain.nii.gz'
ahead_qpd_file = ref_dir+'sub-091_ses-1_acq-wb2_mod-qpd_orient-std_brain.nii.gz'

# main processing steps

if (not os.path.exists(proc_dir)):
    os.makedirs(proc_dir)

print("\n0. MPM clean-up\n")

# skullstrip the MPMs + remove negative values
    
skull = nighres.brain.intensity_based_skullstripping(main_image=mpm_T1w_file, extra_image=mpm_PDw_file,
                    noise_model='exponential', 
                    skip_zero_values=True,
                    iterate=False, dilate_mask=0, topology_lut_dir=None,
                    save_data=True, overwrite=False, output_dir=proc_dir)

brainmask = nighres.io.load_volume(skull['brain_mask']).get_fdata()

mpm_r1strip_file = proc_dir+'sub-091_ses-4_acq-mpmden110_mod-r1_orient-std_brain.nii.gz'
if (not os.path.isfile(mpm_r1strip_file)):
    print("Mask qR1")
    r1 = nighres.io.load_volume(mpm_qr1_file)
    r1strip = nibabel.Nifti1Image(numpy.maximum(0.0,numpy.minimum(2.0,brainmask*r1.get_fdata())), r1.affine, r1.header)
    r1strip = nibabel.as_closest_canonical(r1strip)
    nighres.io.save_volume(mpm_r1strip_file, r1strip)

mpm_r2strip_file = proc_dir+'sub-091_ses-4_acq-mpmden110_mod-r2s_orient-std_brain.nii.gz'
if (not os.path.isfile(mpm_r2strip_file)):
    print("Mask qR2s")
    r2 = nighres.io.load_volume(mpm_qr2s_file)
    r2strip = nibabel.Nifti1Image(numpy.maximum(0.0,numpy.minimum(200.0,brainmask*r2.get_fdata())), r2.affine, r2.header)
    r2strip = nibabel.as_closest_canonical(r2strip)
    nighres.io.save_volume(mpm_r2strip_file, r2strip)

mpm_pdstrip_file = proc_dir+'sub-091_ses-4_acq-mpmden110_mod-pd_orient-std_brain.nii.gz'
if (not os.path.isfile(mpm_pdstrip_file)):
    print("Mask qPD")
    pd = nighres.io.load_volume(mpm_qpd_file)
    pd_data = pd.get_fdata()
    pd_data[pd_data==0] = 200.0
    pd_data[pd_data<0] = 200.0
    pdstrip = nibabel.Nifti1Image(brainmask*pd_data, pd.affine, pd.header)
    pdstrip = nibabel.as_closest_canonical(pdstrip)
    nighres.io.save_volume(mpm_pdstrip_file, pdstrip)


print("\n1. MPM to MP2RAGEME coreg\n")

antsR = nighres.registration.embedded_antspy_multi(
                            source_images=[mpm_r1strip_file,mpm_r2strip_file,mpm_pdstrip_file],
                            target_images=[ahead_qr1_file,ahead_qr2s_file,ahead_qpd_file],
                            run_rigid=True, run_affine=False, run_syn=False,
                            rigid_iterations=10000,
                            affine_iterations=2000,
                            coarse_iterations=180, 
                            medium_iterations=60, fine_iterations=30,
                            cost_function='MutualInformation', 
                            interpolation='Linear',
                            regularization='High',
                            smooth_mask=0.1,
                            ignore_affine=True, 
                            save_data=True, file_name='sub-091_mpmden1102ahead-rigid',
                            output_dir=proc_dir)

antsNL = nighres.registration.embedded_antspy_multi(
                            source_images=[mpm_r1strip_file,mpm_r2strip_file,mpm_pdstrip_file],
                            target_images=[ahead_qr1_file,ahead_qr2s_file,ahead_qpd_file],
                            run_rigid=True, run_affine=True, run_syn=True,
                            rigid_iterations=10000,
                            affine_iterations=2000,
                            coarse_iterations=180, 
                            medium_iterations=60, fine_iterations=30,
                            cost_function='MutualInformation', 
                            interpolation='Linear',
                            regularization='High',
                            smooth_mask=0.1,
                            ignore_affine=True, 
                            save_data=True, file_name='sub-091_mpmden1102ahead-nonlin',
                            output_dir=proc_dir)

antsNL2 = nighres.registration.embedded_antspy_multi(
                            source_images=antsNL['transformed_sources'],
                            target_images=[ahead_qr1_file,ahead_qr2s_file,ahead_qpd_file],
                            run_rigid=True, run_affine=True, run_syn=True,
                            rigid_iterations=10000,
                            affine_iterations=2000,
                            coarse_iterations=180, 
                            medium_iterations=60, fine_iterations=30,
                            cost_function='MutualInformation', 
                            interpolation='Linear',
                            regularization='High',
                            smooth_mask=0.1,
                            ignore_affine=True, 
                            save_data=True, file_name='sub-091_mpmden1102ahead-nonlin2',
                            output_dir=proc_dir)

# map intensities in one step
mpm_r1 = nighres.registration.apply_coordinate_mappings(mpm_r1strip_file,
                            mapping1=antsNL['mapping'],mapping2=antsNL2['mapping'],
                            save_data=True, output_dir=proc_dir)['result']

mpm_r2s = nighres.registration.apply_coordinate_mappings(mpm_r2strip_file,
                            mapping1=antsNL['mapping'],mapping2=antsNL2['mapping'],
                            save_data=True, output_dir=proc_dir)['result']

mpm_pd = nighres.registration.apply_coordinate_mappings(mpm_pdstrip_file,
                            mapping1=antsNL['mapping'],mapping2=antsNL2['mapping'],
                            save_data=True, output_dir=proc_dir)['result']

print("\n2. histogram prior mapping\n")

atlas_dir = '/home/pilou/Projects/Massp-Priors/'
hist = nighres.segmentation.conditional_shape_map_intensities(31, 4, 3,
                          contrast_images=[ahead_qr1_file,ahead_qr2s_file,ahead_qsm_file,ahead_qpd_file], 
                          target_images=[mpm_r1,mpm_r2s,mpm_pd],
                          shape_atlas_probas=atlas_dir+'atlas10_31struct_qmri2fcm_massp-sproba.nii.gz', 
                          shape_atlas_labels=atlas_dir+'atlas10_31struct_qmri2fcm_massp-slabel.nii.gz',
                          intensity_atlas_hist=atlas_dir+'atlas10_31struct_qmri2fcm_massp-chist.nii.gz',
                          skeleton_atlas_probas=atlas_dir+'atlas10_31struct_qmri2fcm_massp-kproba.nii.gz', 
                          skeleton_atlas_labels=atlas_dir+'atlas10_31struct_qmri2fcm_massp-klabel.nii.gz',
                          save_data=True, overwrite=True, output_dir=proc_dir, file_name='atlas10_31struct_mpm_massp-chist.nii.gz')


print("\n3. MASSP for MPM\n")

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
                            save_data=True, file_name='sub-091_ses-4_acq-mpmden110_reg2ahead-step1',
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
                            save_data=True, file_name='sub-091_ses-4_acq-mpmden110_reg2ahead-step2',
                            output_dir=proc_dir)

mapping = nighres.registration.apply_coordinate_mappings(ants1['mapping'],mapping1=ants2['mapping'],
                            save_data=True, overwrite=False, file_name='sub-091_ses-4_acq-mpmden110_reg2ahead-mapping',
                            output_dir=proc_dir)

inverse = nighres.registration.apply_coordinate_mappings(ants2['inverse'],mapping1=ants1['inverse'],
                            save_data=True, overwrite=False, file_name='sub-091_ses-4_acq-mpmden110_reg2ahead-inverse',
                            output_dir=proc_dir)

massp = nighres.parcellation.massp(target_images=[mpm_r1strip_file,mpm_r2strip_file,mpm_pdstrip_file],
                            map_to_target=inverse['result'], 
                            shape_atlas_probas=atlas_dir+'atlas10_31struct_qmri2fcm_massp-sproba.nii.gz', 
                            shape_atlas_labels=atlas_dir+'atlas10_31struct_qmri2fcm_massp-slabel.nii.gz',
                            skeleton_atlas_probas=atlas_dir+'atlas10_31struct_qmri2fcm_massp-kproba.nii.gz', 
                            skeleton_atlas_labels=atlas_dir+'atlas10_31struct_qmri2fcm_massp-klabel.nii.gz',
                            intensity_atlas_hist=proc_dir+'atlas10_31struct_mpm_massp-chist_cspmax-chist.nii.gz',
                            max_iterations=120, max_difference=0.1,
                            intensity_prior=0.25,
                            save_data=True, file_name='sub-091_ses-4_acq-mpmden110',
                            output_dir=proc_dir, overwrite=False)
