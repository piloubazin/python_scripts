import nighres
import numpy
import nibabel
from glob import glob
import os
import subprocess

main_dir='/home/pilou/Projects/Ahead-Database/qMRI-Recomputed/'
in_dir=main_dir
b1_dir=main_dir+'b1maps/'
out_dir=main_dir+'processing/'
final_dir=main_dir

#subjects = range(0,110)
subjects = [98,99,100,101,102,104,105,106,109]

for num in subjects:
    subject = 'sub-'+str(num).zfill(3)
    b1sub = 'Subcortex_'+str(num).zfill(4)+'_'

    mag = sorted(glob(in_dir+subject+'/ses-1/anat/wb/source/'+subject+'*wb*part-mag_mprage.nii.gz'))
    phs = sorted(glob(in_dir+subject+'/ses-1/anat/wb/source/'+subject+'*wb*part-ph_mprage.nii.gz'))
    b1 = sorted(glob(b1_dir+b1sub+'*B1map*.nii.gz'))
                                        
    if (len(mag)==5 and len(phs)==5 and len(b1)==1):         
             
        print("process subject "+subject)

        res_dir = final_dir+subject+'/ses-1/anat/wb/qmri/'    
        r1_final = res_dir+subject+'_ses-1_acq-b1r01_mod-r1hz_orient-std_brain.nii.gz'
        r2_final = res_dir+subject+'_ses-1_acq-b1r01_mod-r2hz_orient-std_brain.nii.gz'
        pd_final = res_dir+subject+'_ses-1_acq-b1r01_mod-qpd_orient-std_brain.nii.gz'
        qsm_final = res_dir+subject+'_ses-1_acq-b1r01_mod-qsm_orient-std_brain.nii.gz'

        print("1. Denoising")
            
        res = nighres.intensity.lcpca_denoising(image_list=mag, phase_list=phs, 
                                          ngb_size=4, stdev_cutoff=1.05,
                                          min_dimension=0, max_dimension=-1,
                                          unwrap=True, process_2d=False,
                                          save_data=True, overwrite=False, output_dir=out_dir)
        
        print("2. B1 mapping")
        
        b1_stack = nighres.io.load_volume(b1[0])
        # extract the 
        b1_img = nibabel.Nifti1Image(b1_stack.get_fdata()[:,:,:,2]/100.0,b1_stack.affine, b1_stack.header)
        b1_file = out_dir+subject+'_b1map-r01.nii.gz'
        nighres.io.save_volume(b1_file,b1_img)
        
        b1w_img = nibabel.Nifti1Image(b1_stack.get_fdata()[:,:,:,1],b1_stack.affine, b1_stack.header)
        b1w_file = out_dir+subject+'_b1w-r01.nii.gz'
        nighres.io.save_volume(b1w_file,b1w_img)
        
        map2wb = nighres.registration.embedded_antsreg(source_image=b1w_file, target_image=res['denoised'][0],
                    run_rigid=True, run_affine=True, run_syn=True,
					ignore_affine=False, ignore_header=False, interpolation="linear", 
                    save_data=True, overwrite=True, output_dir=out_dir)
        
        deform = nighres.registration.apply_coordinate_mappings(b1_file, map2wb['mapping'], 
                                                    interpolation="linear", padding="closest",
                                                    save_data=True, overwrite=True, output_dir=out_dir)
    
        print("3. Quantification")
    
        # whole brain parameters
        qr1 = nighres.intensity.mp2rage_t1_mapping(first_inversion=[res['denoised'][0],                                  
                                        res['denoised'][5]], second_inversion=
                                        [res['denoised'][1],res['denoised'][6]],
                                        correct_B1=True, B1_map=deform['result'], 
                                        inversion_times=[0.670, 3.85],
                                        flip_angles=[4.0, 4.0], inversion_TR=6.72,
                                        excitation_TR=[0.0062, 0.0314],
                                        N_excitations=150, scale_phase=False,
                                        save_data=True,overwrite=False,
                                        output_dir=out_dir, file_name=subject+'_ses1_b1corrected-r01')
        
        qr2 = nighres.intensity.flash_t2s_fitting(image_list=[res['denoised'][1],res['denoised'][2],
                                        res['denoised'][3],res['denoised'][4]],
                                        #te_list=[3.0,11.5,20.0,28.5],
                                        te_list=[0.0030,0.0115,0.0200,0.0285],
                                        save_data=True,overwrite=False,
                                        output_dir=out_dir, file_name=subject+'_ses1_b1corrected')
                                            
        qpd = nighres.intensity.mp2rageme_pd_mapping(first_inversion=[res['denoised'][0],                                  
                                        res['denoised'][5]], second_inversion=
                                        [res['denoised'][1],res['denoised'][6]],
                                        t1map=qr1['t1'], r2smap=qr2['r2s'],
                                        b1map=deform['result'],
                                        echo_times=[0.0030,0.0115,0.0200,0.0285],
                                        inversion_times=[0.670, 3.85],
                                        flip_angles=[4.0, 4.0], inversion_TR=6.72,
                                        excitation_TR=[0.0062, 0.0314],
                                        N_excitations=150,
                                        save_data=True,overwrite=False,
                                        output_dir=out_dir, file_name=subject+'_ses1_b1corrected-r01')
    
        print("3. Skull stripping")
        
        # N4 for the inv2 mag? no difference
        inv2n4_file = os.path.join(out_dir, nighres.utils._fname_4saving(rootfile=res['denoised'][1],
                              suffix='n4'))
        if (not os.path.isfile(inv2n4_file)):
            command = 'N4BiasFieldCorrection -d 3 -i '+res['denoised'][1]+' -o '+inv2n4_file
            print(command)
            try:
                subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
                raise subprocess.CalledProcessError(msg)
        
#        skull = nighres.brain.mp2rage_skullstripping(inv2n4_file, t1_map=qr1['t1'],
#                                        save_data=True, overwrite=False, output_dir=out_dir)
#        skull = nighres.brain.mp2rage_skullstripping(res['denoised'][1], t1_map=qr1['t1'],
#                                        save_data=True, overwrite=False, output_dir=out_dir)
# using s0: not useful
#        skull = nighres.brain.mp2rage_skullstripping(qr2['s0'], t1_map=qr1['t1'],
#                                        save_data=True, overwrite=False, output_dir=out_dir)
        # post processing: tissue mask on inv2
#        fcm = nighres.segmentation.fuzzy_cmeans(skull['inv2_masked'], clusters=5, 
#                    max_difference=0.1, smoothing=0.01, mask_zero=False,
#                    save_data=True, overwrite=False, output_dir=out_dir)

#        fg_file = out_dir+subject+'_ses-1_acq-wb_mod-foreground.nii.gz'
#        if (not os.path.isfile(fg_file)):
#            img = nighres.io.load_volume(fcm['memberships'][0])
#            fg = nibabel.Nifti1Image(1.0-img.get_fdata(), img.affine, img.header)
#            nighres.io.save_volume(fg_file, fg)
        
        skull = nighres.brain.intensity_based_skullstripping(main_image=inv2n4_file, extra_image=qr1['t1'],
                            noise_model='exponential', 
                            skip_zero_values=True,
                            iterate=False, dilate_mask=1, topology_lut_dir=None,
                            save_data=True, overwrite=False, output_dir=out_dir)

#        mask = nighres.io.load_volume(skull['brain_mask']).get_fdata()
#        fg = nighres.io.load_volume(fcm['classification']).get_fdata()>1
#        brainmask = fg*mask
        
#        finalmask_file = out_dir+subject+'_ses-1_acq-wb_mod-brainmask.nii.gz'
#        if (not os.path.isfile(finalmask_file)):
#            img = nighres.io.load_volume(skull['brain_mask'])
#            finalmask = nibabel.Nifti1Image(brainmask, img.affine, img.header)
#            finalmask = nibabel.as_closest_canonical(finalmask)
#            nighres.io.save_volume(finalmask_file, finalmask)

        # stripping other contrasts? R1, R2*, PD
        # and rescale to Hz for R1, R2*
        brainmask = nighres.io.load_volume(skull['brain_mask']).get_fdata()
        
        r1strip_file = out_dir+subject+'_ses-1_acq-b1r01_mod_r1hz_orient-std_brain.nii.gz'
        if (not os.path.isfile(r1strip_file) and not os.path.isfile(r1_final)):
            print("Mask qR1")
            r1 = nighres.io.load_volume(qr1['r1'])
            r1strip = nibabel.Nifti1Image(numpy.minimum(3.0,brainmask*r1.get_fdata()), r1.affine, r1.header)
            r1strip = nibabel.as_closest_canonical(r1strip)
            nighres.io.save_volume(r1strip_file, r1strip)
    
        r2strip_file = out_dir+subject+'_ses-1_acq-b1r01_mod_r2hz_orient-std_brain.nii.gz'
        if (not os.path.isfile(r2strip_file) and not os.path.isfile(r2_final)):
            print("Mask qR2*")
            r2s = nighres.io.load_volume(qr2['r2s'])
            r2strip = nibabel.Nifti1Image(numpy.minimum(200.0,brainmask*r2s.get_fdata()), r2s.affine, r2s.header)
            r2strip = nibabel.as_closest_canonical(r2strip)
            nighres.io.save_volume(r2strip_file, r2strip)
    
        # call N4
        pdn4_file = os.path.join(out_dir, nighres.utils._fname_4saving(rootfile=qpd['pd'],
                              suffix='n4'))
        if (not os.path.isfile(pdn4_file)):
            command = 'N4BiasFieldCorrection -d 3 -i '+qpd['pd']+'  -x '+skull['brain_mask']+' -o '+pdn4_file
            print(command)
            try:
                subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
                raise subprocess.CalledProcessError(msg)
        
        pdstrip_file = out_dir+subject+'_ses-1_acq-b1r01_mod_qpd_orient-std_brain.nii.gz'
        if (not os.path.isfile(pdstrip_file) and not os.path.isfile(pd_final)):
            print("Mask qPD")
            pd = nighres.io.load_volume(pdn4_file)
            pddata = pd.get_fdata()
            pdmean = numpy.mean(pddata[pddata>0])
            pdstrip = nibabel.Nifti1Image(numpy.minimum(8.0,brainmask*pddata/pdmean), pd.affine, pd.header)
            pdstrip = nibabel.as_closest_canonical(pdstrip)
            nighres.io.save_volume(pdstrip_file, pdstrip)
    
        print("4. QSM")
            
        magnitude = res['denoised'][1]
        phase2 = res['denoised'][7]
        phase3 = res['denoised'][8]
        phase4 = res['denoised'][9]
        phases = [phase2,phase3,phase4]
        tes = [11.5,20.0,28.5]
        
        for idx,phs in enumerate(phases):
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
                
                # here use the brain mask from earlier, but copy the header from phase
                # (the N4 step may have mangled the mask header...)
                #mask = skull['brain_mask']
                qsm_mask_file = os.path.join(out_dir, nighres.utils._fname_4saving(rootfile=phs,
                              suffix='tgv-mask'))
                phs_img = nighres.io.load_volume(phs)
                mask_img = nibabel.Nifti1Image(brainmask, phs_img.affine, phs_img.header)
                nighres.io.save_volume(qsm_mask_file, mask_img)
                mask = qsm_mask_file
                
                # rescale if needed? not after denoising
                
                # run tgv_qsm on each echo separately at first
                recon = 'tgv_qsm -p '+phs
                #recon = recon + ' -m '+mask['mask']
                recon = recon + ' -m '+mask
                recon = recon + ' -e '+str(5)
                recon = recon + ' -i '+str(1000)
                recon = recon + ' -f '+str(7)
                recon = recon + ' -t '+str(0.001*tes[idx])
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
        #qsm_file = out_dir+subject+'_ses-1_acq-wb_mod_qsm_orient-std_brain.nii.gz'
        qsm_file = os.path.join(out_dir, nighres.utils._fname_4saving(rootfile=phases[0],
                              suffix='tgv-qsm_med-img'))
        
        rng_file = os.path.join(out_dir, nighres.utils._fname_4saving(rootfile=phases[0],
                              suffix='tgv-qsm_rng-img'))
        
        if ((not os.path.isfile(qsm_file) and not os.path.isfile(qsm_final) ) \
            or not os.path.isfile(rng_file)):
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
            med = numpy.percentile(data, 50, axis=3)
            rng = data.max(axis=3)-data.min(axis=3)
            
            img = nighres.io.load_volume(qsm2)
            img = nibabel.Nifti1Image(numpy.maximum(-0.3,numpy.minimum(0.5,brainmask*med)),img.affine,img.header)
            img = nibabel.as_closest_canonical(img)
            nighres.io.save_volume(qsm_file, img)
    
            img = nibabel.Nifti1Image(brainmask*rng,img.affine,img.header)
            nighres.io.save_volume(rng_file, img)
    
        print("5. Copy final images to BIDS folder")
            
        if (not os.path.exists(res_dir)):
            os.makedirs(res_dir)
            
        if (not os.path.exists(r1_final)):
            os.rename(r1strip_file,r1_final)
            
        if (not os.path.exists(r2_final)):            
            os.rename(r2strip_file,r2_final)
            
        if (not os.path.exists(pd_final)):
            os.rename(pdstrip_file,pd_final)
                        
        if (not os.path.exists(qsm_final)):
            os.rename(qsm_file,qsm_final)     
