import nighres
import numpy
import nibabel
from glob import glob
import os
import subprocess

main_dir='/home/pilou/Projects/Ahead-Database/Vascular-Modeling/'
in_dir=main_dir+'fatnav-corrected/'
out_dir=main_dir+'processing/'
final_dir=main_dir+'vessel-maps/'
subjects = ['sub-004','sub-006']

for subject in subjects:

    print("process subject "+subject)
    
    print("1. Denoising")
    
    mag = sorted(glob(in_dir+subject+'_*_m_corr.nii.gz'))
    phs = sorted(glob(in_dir+subject+'_*_ph_corr.nii.gz'))
    
    res = nighres.intensity.lcpca_denoising(image_list=mag, phase_list=phs, 
                                      ngb_size=4, stdev_cutoff=1.05,
                                      min_dimension=0, max_dimension=-1,
                                      unwrap=True, process_2d=False,
                                      save_data=True, overwrite=False, output_dir=out_dir)
    
    print("2. Quantification")

# slab parameters
#    qr1 = nighres.intensity.mp2rage_t1_mapping(first_inversion=[res['denoised'][0],
#                                    res['denoised'][5]], second_inversion=
#                                    [res['denoised'][1],res['denoised'][6]],
#                                    inversion_times=[0.670, 3.7377],
#                                    flip_angles=[7.0, 8.0], inversion_TR=8.330,
#                                    excitation_TR=[0.008, 0.032],
#                                    N_excitations=150,
#                                    save_data=True,overwrite=False,
#                                    output_dir=out_dir)
    
#    qr2 = nighres.intensity.flash_t2s_fitting(image_list=[res['denoised'][1],res['denoised'][2],
#                                    res['denoised'][3],res['denoised'][4]],
#                                    te_list=[4.6,12.6,20.6,28.6],
#                                    save_data=True,overwrite=False,
#                                    output_dir=out_dir)

    # whole brain parameters
    qr1 = nighres.intensity.mp2rage_t1_mapping(first_inversion=[res['denoised'][0],                                  
                                    res['denoised'][5]], second_inversion=
                                    [res['denoised'][1],res['denoised'][6]],
                                    inversion_times=[0.670, 3.85],
                                    flip_angles=[7.0, 6.0], inversion_TR=6.72,
                                    excitation_TR=[0.0062, 0.0314],
                                    N_excitations=150,
                                    save_data=True,overwrite=False,
                                    output_dir=out_dir)
    
    qr2 = nighres.intensity.flash_t2s_fitting(image_list=[res['denoised'][1],res['denoised'][2],
                                    res['denoised'][3],res['denoised'][4]],
                                    te_list=[3.0,11.5,20.0,28.5],
                                    save_data=True,overwrite=False,
                                    output_dir=out_dir)
                                        
    print("3. Skull stripping")
    
    skull = nighres.brain.mp2rage_skullstripping(res['denoised'][1], t1_map=qr1['t1'],
                                    save_data=True, overwrite=False, output_dir=out_dir)
    
    
    print("4. QSM")
        
    magnitude = res['denoised'][1]
    phase2 = res['denoised'][7]
    phase3 = res['denoised'][8]
    phase4 = res['denoised'][9]
    phases = [phase2,phase3,phase4]
    
    for phs in phases:
        # check if already computed
        qsm_file = os.path.join(out_dir, nighres.utils._fname_4saving(rootfile=phs,
                          suffix='tgv-qsm'))
    
        if (not os.path.isfile(qsm_file)):
            # estimate background from magnitude
            #mask = nighres.intensity.background_estimation(image = magnitude, 
            #              distribution='exponential', ratio=1e-3,
            #              skip_zero=True, iterate=False, dilate=-1,
            #              threshold=0.5, overwrite=True,
            #              save_data=True, output_dir=out_dir)
            
            # here use the brain mask from earlier
            mask = skull['brain_mask']
            
            # rescale if needed? not after denoising
            
            # run tgv_qsm on each echo separately at first
            recon = 'tgv_qsm -p '+phs
            #recon = recon + ' -m '+mask['mask']
            recon = recon + ' -m '+mask
            recon = recon + ' -i '+str(1000)
            recon = recon + ' -f '+str(7)
            recon = recon + ' -t '+str(0.0286)
            recon = recon + ' --alpha '+str(0.0015)+' '+str(0.0005)
            recon = recon + ' -o '+'_reconQSM'
            recon = recon + ' -v'
            
            print(recon)
            try:
                subprocess.check_output(recon, shell=True, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                msg = 'execution failed (error code '+e.returncode+')\n Output: '+e.output
                raise subprocess.CalledProcessError(msg)
            
            # move the results in proper folder
            results = sorted(glob(out_dir+subject+'*_reconQSM*.nii.gz'))
            for result in results:
                os.rename(result, qsm_file)
    
    # recombine the three echos
    qsm_file = os.path.join(out_dir, nighres.utils._fname_4saving(rootfile=phases[0],
                          suffix='tgv-qsm_med-img'))
    
    if (not os.path.isfile(qsm_file)):
        qsm2 = os.path.join(out_dir, nighres.utils._fname_4saving(rootfile=phase2,
                          suffix='tgv-qsm'))
        qsm3 = os.path.join(out_dir, nighres.utils._fname_4saving(rootfile=phase3,
                          suffix='tgv-qsm'))
        qsm4 = os.path.join(out_dir, nighres.utils._fname_4saving(rootfile=phase4,
                          suffix='tgv-qsm'))
        
        data2 = nighres.io.load_volume(qsm2).get_data()
        data3 = nighres.io.load_volume(qsm3).get_data()
        data4 = nighres.io.load_volume(qsm4).get_data()
        
        data = numpy.stack((data2,data3,data4),axis=3)
        data = numpy.percentile(data, 50, axis=3)
        
        img = nighres.io.load_volume(qsm2)
        img = nibabel.Nifti1Image(data,img.affine,img.header)
        nighres.io.save_volume(qsm_file, img)

    # 5. Post-processing: subcortex, cortex, vasculature segmentation (...): see Pilou

