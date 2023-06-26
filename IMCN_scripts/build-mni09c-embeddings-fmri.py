import nibabel
import nighres
import glob
import numpy
import os
import sys
            
in_dir = '/home/pilou/Projects/Ahead-Database/Coordinate-Systems/'
out_dir = in_dir+'flatmaps_fmri/'

atlas_dir = '/home/pilou/Datasets/mni_icbm152_nlin_asym_09c/'
mni = atlas_dir+'mni_icbm152_t1_tal_nlin_asym_09c.nii.gz'
wm = atlas_dir+'mni_icbm152_wm_tal_nlin_asym_09c.nii.gz'
gm = atlas_dir+'mni_icbm152_gm_tal_nlin_asym_09c.nii.gz'

massp_dir = '/home/pilou/Projects/Ahead-Database/Automated-Parcellation/atlas_maps/qmri2fcm2021/massp2021_subcortex_lifespan/'
massp = massp_dir+'massp2021-parcellation_decade-18to80.nii.gz'

massp_rm = [15,16,17,18,21,22,25]

template_dir = '/home/pilou/Projects/Ahead-Database/MNI05-Coregistration/mni_icbm152_nlin_asym_09b/'
template = template_dir+'mni_icbm152_t1_tal_nlin_asym_09b_hires.nii.gz'

func_dir = '/home/pilou/Projects/Ahead-Database/Coordinate-Systems/func_maps/'
func_map1 = func_dir+'model-0a-z_fwhm-1p5_subjectlevelcontrast-0_grouplevelcontrast-1_flame-12_desc-zstat.nii.gz'
func_map2 = func_dir+'model-0a-z_fwhm-1p5_subjectlevelcontrast-5_grouplevelcontrast-1_flame-12_desc-zstat.nii.gz'

# coreg MASSP labels to MNI2009c
ants = nighres.registration.embedded_antspy(
                                    source_image=template,
                                    target_image=mni,
                                    run_rigid=True, run_affine=False, run_syn=False,
                                    rigid_iterations=10000,
                                    affine_iterations=2000,
                                    coarse_iterations=180, 
                                    medium_iterations=60, fine_iterations=30,
                                    cost_function='MutualInformation', 
                                    interpolation='Linear',
                                    regularization='High',
                                    smooth_mask=0.0,
                                    ignore_affine=False, 
                                    save_data=True, file_name='mni09b2mni09c',
                                    output_dir=out_dir)

massp_img = nighres.io.load_volume(massp)

massp09c = nighres.registration.apply_coordinate_mappings(massp_img,mapping1=ants['mapping'],
                                    interpolation='nearest',
                                    save_data=True, overwrite=False, file_name='massp2mni09c',
                                    output_dir=out_dir)

# build Laplacian map

wmgm = out_dir+'mni_icbm152_wmgm_tal_nlin_asym_09c.nii.gz'
if not os.path.exists(wmgm):
    wm_img = nighres.io.load_volume(wm)
    gm_img = nighres.io.load_volume(gm)

    wmgm_img = nibabel.Nifti1Image(wm_img.get_fdata()+gm_img.get_fdata(),wm_img.affine,wm_img.header)
    nighres.io.save_volume(wmgm, wmgm_img)
    
thr = nighres.io.load_volume(wmgm)
thr = nibabel.Nifti1Image(thr.get_fdata()>0.95,thr.affine,thr.header)
                
massp_img = nighres.io.load_volume(massp09c['result']) 

massp_data = massp_img.get_fdata()
for lb in massp_rm:
    massp_data[massp_data==lb] = 0    
massp_img = nibabel.Nifti1Image(massp_data,massp_img.affine,massp_img.header)
nighres.io.save_volume(out_dir+'massp2mni09c_gm.nii.gz', massp_img)
                
                
redo = False                
                
shape = nighres.shape.spectral_embedding(thr,msize=800,contrasts=[massp_img],
                                    scaling=0.0,factor=5.0,alpha=0.0,
                                    bg="boundary",step=0.01, 
                                    save_data=True,overwrite=redo,output_dir=out_dir,
                                    file_name='mni09c_wmgm_th095_massp_f5_bound')
                
                

# project fMRI results

maps = nighres.shape.spectral_flatmap(thr,shape['result'],size=512,dims=2,combined=False,
                                    save_data=True,output_dir=out_dir,overwrite=redo,
                                    file_name='mni09c_wmgm_th095_massp_f5_bound_ax512')
    
maps = nighres.shape.spectral_flatmap(thr,shape['result'],contrast_image=massp_img,
                                    size=512,dims=2,offset=0,combined=False,contrast_mode="min-bound",
                                    save_data=True,output_dir=out_dir,overwrite=True,
                                    file_name='mni09c_wmgm_th095_massp_f5_bound_massp512')

maps = nighres.shape.spectral_flatmap(thr,shape['result'],contrast_image=func_map1,
                                    size=512,dims=2,offset=0,combined=False,
                                    save_data=True,output_dir=out_dir,overwrite=redo,
                                    file_name='mni09c_wmgm_th095_massp_f5_bound_fmri0512')

maps = nighres.shape.spectral_flatmap(thr,shape['result'],contrast_image=func_map2,
                                    size=512,dims=2,offset=0,combined=False,
                                    save_data=True,output_dir=out_dir,overwrite=redo,
                                    file_name='mni09c_wmgm_th095_massp_f5_bound_fmri5512')


maps = nighres.shape.spectral_flatmap(thr,shape['result'],size=512,dims=2,offset=1,combined=False,
                                    save_data=True,output_dir=out_dir,overwrite=redo,
                                    file_name='mni09c_wmgm_th095_massp_f5_bound_cr512')
    
maps = nighres.shape.spectral_flatmap(thr,shape['result'],contrast_image=massp_img,
                                    size=512,dims=2,offset=1,combined=False,contrast_mode="min-bound",
                                    save_data=True,output_dir=out_dir,overwrite=True,
                                    file_name='mni09c_wmgm_th095_massp_f5_bound_masspcr512')

maps = nighres.shape.spectral_flatmap(thr,shape['result'],contrast_image=func_map1,
                                    size=512,dims=2,offset=1,combined=False,
                                    save_data=True,output_dir=out_dir,overwrite=redo,
                                    file_name='mni09c_wmgm_th095_massp_f5_bound_fmri0cr512')

maps = nighres.shape.spectral_flatmap(thr,shape['result'],contrast_image=func_map2,
                                    size=512,dims=2,offset=1,combined=False,
                                    save_data=True,output_dir=out_dir,overwrite=redo,
                                    file_name='mni09c_wmgm_th095_massp_f5_bound_fmri5cr512')


maps = nighres.shape.spectral_flatmap(thr,shape['result'],size=512,dims=2,offset=2,combined=False,
                                    save_data=True,output_dir=out_dir,overwrite=redo,
                                    file_name='mni09c_wmgm_th095_massp_f5_bound_sg512')
    
maps = nighres.shape.spectral_flatmap(thr,shape['result'],contrast_image=massp_img,
                                    size=512,dims=2,offset=2,combined=False,contrast_mode="min-bound",
                                    save_data=True,output_dir=out_dir,overwrite=True,
                                    file_name='mni09c_wmgm_th095_massp_f5_bound_masspsg512')

maps = nighres.shape.spectral_flatmap(thr,shape['result'],contrast_image=func_map1,
                                    size=512,dims=2,offset=2,combined=False,
                                    save_data=True,output_dir=out_dir,overwrite=redo,
                                    file_name='mni09c_wmgm_th095_massp_f5_bound_fmri0sg512')

maps = nighres.shape.spectral_flatmap(thr,shape['result'],contrast_image=func_map2,
                                    size=512,dims=2,offset=2,combined=False,
                                    save_data=True,output_dir=out_dir,overwrite=redo,
                                    file_name='mni09c_wmgm_th095_massp_f5_bound_fmri5sg512')


