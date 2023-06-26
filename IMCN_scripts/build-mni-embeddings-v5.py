import nibabel
import nighres
import glob
import numpy
import os
import sys
            
in_dir = '/home/pilou/Projects/Ahead-Database/Coordinate-Systems/'
out_dir = in_dir+'flatmaps_allbg/'

mni1 = in_dir+'mni_icbm152_gm_wm_tal_nlin_asym_09a.nii.gz'
mni2 = in_dir+'mni_icbm152_gm_wm_tal_nlin_asym_09a_sub-2.nii.gz'
mni4 = in_dir+'mni_icbm152_gm_wm_tal_nlin_asym_09a_sub-4.nii.gz'

tis1 = in_dir+'mni_icbm152_tissues_tal_nlin_asym_09a.nii.gz'
tis2 = in_dir+'mni_icbm152_tissues_tal_nlin_asym_09a_sub-2.nii.gz'
tis4 = in_dir+'mni_icbm152_tissues_tal_nlin_asym_09a_sub-4.nii.gz'

wm1 = in_dir+'mni_icbm152_wm_tal_nlin_asym_09a.nii.gz'
wm2 = in_dir+'mni_icbm152_wm_tal_nlin_asym_09a_sub-2.nii.gz'
wm4 = in_dir+'mni_icbm152_wm_tal_nlin_asym_09a_sub-4.nii.gz'

gm1 = in_dir+'mni_icbm152_gm_tal_nlin_asym_09a.nii.gz'
gm2 = in_dir+'mni_icbm152_gm_tal_nlin_asym_09a_sub-2.nii.gz'
gm4 = in_dir+'mni_icbm152_gm_tal_nlin_asym_09a_sub-4.nii.gz'

massp1 = in_dir+'massp2mni2_def-img_gm.nii.gz'
massp2 = in_dir+'massp2mni1_def-img_gm.nii.gz'
massp4 = in_dir+'massp2mni0_def-img_gm.nii.gz'

mni1massp = in_dir+'mni_icbm152_gm_wm_tal_nlin_asym_09a_wmassp.nii.gz'
massp1full = in_dir+'massp2mni2_def-img.nii.gz'

ding1 = in_dir+'annotation_full_1mm.nii.gz'
ding2 = in_dir+'annotation_full_2mm.nii.gz'
ding4 = in_dir+'annotation_full_4mm.nii.gz'


mni = [mni4,mni2,mni1]
tis = [tis4,tis2,tis1]
gm = [gm4,gm2,gm1]
wm = [wm4,wm2,wm1]
massp = [massp4,massp2,massp1]
ding = [ding4,ding2,ding1]

#debug
mni = [mni4,mni2,mni1]


# weighting by MASSP
for idx,img in enumerate(mni):
    for factor in [2.0,3.0,4.0]:
        for alpha in [0.0]:
            for scaling in [-100.0,0.0,100.0]:
                redo=False
                print("process "+img)
                print(" with contrast "+massp[idx])
                
                thr = nighres.io.load_volume(img)
                thr = nibabel.Nifti1Image(thr.get_fdata()>0.95,thr.affine,thr.header)
            
                wgt = massp[idx];
            
                shape = nighres.shape.spectral_embedding(thr,msize=800,contrasts=[wgt],scaling=scaling,factor=factor,alpha=alpha,bg="boundary",
                                             step=0.01,
                                             save_data=True,overwrite=redo,output_dir=out_dir,
                                             file_name='mni_wmgm_th095_masspAcontrast_sc-'+str(idx)+'_f-'+str(int(factor))+'_a-'+str(int(10*alpha))+'_d-'+str(int(10*scaling)))
            
    #            shape = nighres.shape.spectral_tsne(thr,shape['result'],contrasts=[wgt],scaling=1.0,factor=factor,alpha=alpha,bg="boundary",
    #                                         step=200.0,momentum=0.0,relaxation=0.0,iterations=200,
    #                                         save_data=True,overwrite=redo,output_dir=out_dir,
    #                                         file_name='mni_wmgm_th095_masspAcontrast_sc-'+str(idx)+'_f-'+str(int(factor))+'_a-'+str(int(10*alpha))+'_tsne.nii.gz')
            
                maps = nighres.shape.spectral_flatmap(thr,shape['result'],size=512,dims=2,combined=False,
                                            save_data=True,output_dir=out_dir,overwrite=redo,
                                            file_name='mni_wmgm_th095_masspAcontrast_sc-'+str(idx)+'_f-'+str(int(factor))+'_a-'+str(int(10*alpha))+'_d-'+str(int(10*scaling))+'_ax.nii.gz')
    
                maps = nighres.shape.spectral_flatmap(thr,shape['result'],
                                            contrast_image=massp[idx],
                                            size=512,dims=2,offset=0,combined=False,
                                            save_data=True,output_dir=out_dir,overwrite=redo,
                                            file_name='mni_wmgm_th095_masspAcontrast_sc-'+str(idx)+'_f-'+str(int(factor))+'_a-'+str(int(10*alpha))+'_d-'+str(int(10*scaling))+'_massp.nii.gz')
