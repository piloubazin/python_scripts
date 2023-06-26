import nighres
import numpy
import nibabel
from glob import glob
import os
import subprocess
import math
import scipy.ndimage
import ants

# hopefully we just have to change the subject id every time
#subjects = ['sub-010','sub-011','sub-012','sub-013','sub-014','sub-018','sub-019']
#subjects = ['sub-038','sub-030','sub-033','sub-034']
#subjects = ['sub-031','sub-032','sub-035','sub-039','sub-040','sub-041','sub-042']
subjects = ['sub-043','sub-044']
session = 'ses-anatomical'

for subject in subjects:

    # in case there was motion between sessions
    coalign = True
    
    # keywords identifying the three different sequences.
    # they shouldn't change except for typos during scanning
    #gre_name = 'ke_gre_aspire'
    #gre_mag_number = '17'
    #gre_phs_number = '16'
    gre_mag = 'Aspire_M_ke_gre_aspire'
    gre_phs = 'Aspire_P_ke_gre_aspire'
    
    mp2rage_name = 't1_mp2rage'
    mp2rage_inv1 = 't1_mp2rage_sag_0p75_INV1'
    mp2rage_inv2 = 't1_mp2rage_sag_0p75_INV2'
    mp2rage_uni = 't1_mp2rage_sag_0p75_UNI'
    
    b1map_name = 'b1map'
    
    main_dir='/home/Public/trondheim/'
    
    in_file=main_dir+'sourcedata/zipdata/'+subject+'/'+session+'/'+subject+'_'+session+'_data.zip'
    proc_dir=main_dir+'processing/'+subject+'/'
    res_dir=main_dir+'derivatives/nighres/'+subject+'/'+session+'/qmri/'
    
    # output names
    final_qr1_file = res_dir+subject+'_'+session+'_acq-mp2rage_mod-r1hz_orient-std_brain.nii.gz'
    final_qr2s_file = res_dir+subject+'_'+session+'_acq-gre_mod-r2hz_orient-std_brain.nii.gz'
    final_qsm_file = res_dir+subject+'_'+session+'_acq-gre_mod-qsm_orient-std_brain.nii.gz'
    final_qpd_file = res_dir+subject+'_'+session+'_acq-mp2rage_mod-qpd_orient-std_brain.nii.gz'
    
    final_qr1nob1_file = res_dir+subject+'_'+session+'_acq-mp2rage_mod-r1hznob1_orient-std_brain.nii.gz'
    final_qpdnob1_file = res_dir+subject+'_'+session+'_acq-mp2rage_mod-qpdnob1_orient-std_brain.nii.gz'
    
    # here we look for data, assuming names are consistent
    
    print("\n0. Uncompressing and Nifti conversion\n")
    
    if (not os.path.exists(proc_dir)):
        os.makedirs(proc_dir)
        
    if (not os.path.exists(proc_dir+subject+'_'+session+'_data/SCANS/')):
        command = 'unzip '+in_file+' -d '+proc_dir
        print(command)
        try:
            subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
            raise subprocess.CalledProcessError(msg)
    
    if (len(glob(proc_dir+'*.nii.gz'))==0):    
        #command = 'dcm2niix -z y -o '+proc_dir+' -f %p_%s_e%e '+proc_dir+subject+'_'+session+'_data/SCANS/'
        command = 'dcm2niix -z y -o '+proc_dir+' -f %d_%s_e%e '+proc_dir+subject+'_'+session+'_data/SCANS/'
        print(command)
        try:
            subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
            raise subprocess.CalledProcessError(msg)
    
    
    
    print("\n1. Looking for data\n")
    
    # multi-echo GRE (aspire, first image is not used)
    #gre_mag_e1 = sorted(glob(proc_dir+'*'+gre_name+'*'+gre_mag_number+'_e1.nii*'))
    #gre_mag_e2 = sorted(glob(proc_dir+'*'+gre_name+'*'+gre_mag_number+'_e2.nii*'))
    #gre_mag_e3 = sorted(glob(proc_dir+'*'+gre_name+'*'+gre_mag_number+'_e3.nii*'))
    #gre_mag_e4 = sorted(glob(proc_dir+'*'+gre_name+'*'+gre_mag_number+'_e4.nii*'))
    
    #gre_phs_e1 = glob(proc_dir+'*'+gre_name+'*'+gre_phs_number+'_e1_ph.nii*')
    #gre_phs_e2 = glob(proc_dir+'*'+gre_name+'*'+gre_phs_number+'_e2_ph.nii*')
    #gre_phs_e3 = glob(proc_dir+'*'+gre_name+'*'+gre_phs_number+'_e3_ph.nii*')
    #gre_phs_e4 = glob(proc_dir+'*'+gre_name+'*'+gre_phs_number+'_e4_ph.nii*')
    
    gre_mag_e1 = sorted(glob(proc_dir+'*'+gre_mag+'*_e1.nii*'))
    gre_mag_e2 = sorted(glob(proc_dir+'*'+gre_mag+'*_e2.nii*'))
    gre_mag_e3 = sorted(glob(proc_dir+'*'+gre_mag+'*_e3.nii*'))
    gre_mag_e4 = sorted(glob(proc_dir+'*'+gre_mag+'*_e4.nii*'))
    
    gre_phs_e1 = glob(proc_dir+'*'+gre_phs+'*_e1_ph.nii*')
    gre_phs_e2 = glob(proc_dir+'*'+gre_phs+'*_e2_ph.nii*')
    gre_phs_e3 = glob(proc_dir+'*'+gre_phs+'*_e3_ph.nii*')
    gre_phs_e4 = glob(proc_dir+'*'+gre_phs+'*_e4_ph.nii*')
    
    if (len(gre_mag_e1)>0 and len(gre_mag_e2)>0 and len(gre_mag_e3)>0 and len(gre_mag_e4)>0 \
        and len(gre_phs_e1)>0 and len(gre_phs_e2)>0 and len(gre_phs_e3)>0 and len(gre_phs_e4)>0) :
    #    print("Aspire GRE found. Using images "+gre_phs_number+" and "+gre_mag_number+" as phase and magnitude\n")
        print("Aspire GRE found. Using images "+gre_phs+" and "+gre_mag+" as phase and magnitude\n")
        
        gre_mag_e1 = gre_mag_e1[0]
        gre_mag_e2 = gre_mag_e2[0]
        gre_mag_e3 = gre_mag_e3[0]
        gre_mag_e4 = gre_mag_e4[0]
    
        gre_phs_e1 = gre_phs_e1[0]
        gre_phs_e2 = gre_phs_e2[0]
        gre_phs_e3 = gre_phs_e3[0]
        gre_phs_e4 = gre_phs_e4[0]
    
    # MP2RAGE
    uni = sorted(glob(proc_dir+mp2rage_uni+'*.nii*'))
    inv1 = sorted(glob(proc_dir+mp2rage_inv1+'*.nii*'))
    inv2 = sorted(glob(proc_dir+mp2rage_inv2+'*.nii*'))
    
    if (len(uni)>0 and len(inv1)>0 and len(inv2)>0) :
        print("Mp2rage found.")
        uni = uni[0]
        inv1 = inv1[0]
        inv2 = inv2[0]
    else:
        mp2rage = sorted(glob(proc_dir+'*'+mp2rage_name+'*.nii*'))
        
        if len(mp2rage)>2:
            print("MP2RAGE found. Checking intensity histograms to set INV1 INV2 and UNI\n")
            
            max0 = numpy.max(nighres.io.load_volume(mp2rage[0]).get_fdata())
            max1 = numpy.max(nighres.io.load_volume(mp2rage[1]).get_fdata())
            max2 = numpy.max(nighres.io.load_volume(mp2rage[2]).get_fdata())
            
            if (max0>max1 and max0>max2 and max2>max1):
                uni = mp2rage[0]
                inv1 = mp2rage[1]
                inv2 = mp2rage[2]
            elif (max0>max1 and max0>max2 and max1>max2):
                uni = mp2rage[0]
                inv1 = mp2rage[2]
                inv2 = mp2rage[1]
            elif (max1>max0 and max1>max2 and max2>max0):
                uni = mp2rage[1]
                inv1 = mp2rage[0]
                inv2 = mp2rage[2]
            elif (max1>max0 and max1>max2 and max0>max2):
                uni = mp2rage[1]
                inv1 = mp2rage[2]
                inv2 = mp2rage[0]
            elif (max2>max0 and max2>max1 and max1>max0):
                uni = mp2rage[2]
                inv1 = mp2rage[0]
                inv2 = mp2rage[1]
            elif (max2>max0 and max2>max1 and max0>max1):
                uni = mp2rage[2]
                inv1 = mp2rage[1]
                inv2 = mp2rage[0]
     
    # B1 map (optional)
    b1map = sorted(glob(proc_dir+'*'+b1map_name+'*.nii*'))
    
    b1img = None
    b1ratio = None
    if len(b1map)>1:
        print("B1 map found. Using first as anatomy, second as ratio\n")
    
        b1img = b1map[0]
        b1ratio = b1map[1]
    
    # if R2* not acquired with same geometry as MP2RAGE
    inv2_img = nighres.io.load_volume(inv2)
    gre1_img = nighres.io.load_volume(gre_mag_e1)
    
    if ( (inv2_img.shape is not gre1_img.shape) or coalign):
        print ("Warning: MP2RAGE and GRE have different dimensions (and/or may have moved), using scanner info for alignment")
        print ("MP2RAGE is used as anatomical basis, but GRE-based mapping is done in native space (interpolation issues)\n")
        
        mapping = nighres.registration.embedded_antspy(gre_mag_e1, inv2,
                    run_rigid=coalign,
                    run_affine=False,
                    run_syn=False,
                    interpolation='Linear',
                    ignore_affine=False, ignore_header=False,
                    save_data=True, overwrite=False, output_dir=proc_dir)
        
    print("\n2. GRE Denoising\n")
    
    # run the denoising on gre only
    mag = [gre_mag_e1,gre_mag_e2,gre_mag_e3,gre_mag_e4]
    phs = [gre_phs_e1,gre_phs_e2,gre_phs_e3,gre_phs_e4]
    
    res = nighres.intensity.lcpca_denoising(image_list=mag, phase_list=phs, 
                                  ngb_size=4, stdev_cutoff=1.05,
                                  min_dimension=0, max_dimension=-1,
                                  unwrap=True, rescale_phs=True, process_2d=False,
                                  save_data=True, overwrite=False, output_dir=proc_dir)
    
    gre_mag_e1 = res['denoised'][0]
    gre_mag_e2 = res['denoised'][1]
    gre_mag_e3 = res['denoised'][2]
    gre_mag_e4 = res['denoised'][3]
    
    gre_phs_e1 = res['denoised'][4]
    gre_phs_e2 = res['denoised'][5]
    gre_phs_e3 = res['denoised'][6]
    gre_phs_e4 = res['denoised'][7]
    
    
        
    print("\n3. T1, T2*, PD mapping\n")
    
    # basic T2* mapping
    qr2 = nighres.intensity.flash_t2s_fitting(image_list=[gre_mag_e1,gre_mag_e2,gre_mag_e3,gre_mag_e4],
                                    te_list=[0.00251,0.00722,0.01444,0.02323],
                                    save_data=True,overwrite=False,
                                    output_dir=proc_dir)              

    if ( (inv2_img.shape is not gre1_img.shape) or coalign):
        # transform the T2 map to T1 space after processing
        qr2map = nighres.registration.apply_coordinate_mappings(qr2['r2s'], mapping['mapping'],
                        interpolation="nearest", padding="zero",
                        save_data=True, overwrite=False, output_dir=proc_dir)['result']
    else:
        qr2map = qr2['r2s']    

    # T1 and PD mapping
    if b1ratio is not None:
        b1map = nighres.registration.embedded_antspy(b1ratio, uni,
                    run_rigid=False,
                    run_affine=False,
                    run_syn=False,
                    interpolation='Linear',
                    ignore_affine=False, ignore_header=False,
                    save_data=True, overwrite=False, output_dir=proc_dir)
        
        qr1 = nighres.intensity.mp2rage_t1_from_uni(uniform_image=uni,
                                        inversion_times=[0.840, 2.370],
                                        flip_angles=[5.0, 6.0], inversion_TR=4.3,
                                        excitation_TR=[0.0072, 0.0072],
                                        N_excitations=180,
                                        correct_B1=True, B1_map=b1map['transformed_source'], 
                                        B1_scale=1000.0,
                                        save_data=True,overwrite=False,
                                        output_dir=proc_dir)

        qpd = nighres.intensity.mp2rageme_pd_mapping(first_inversion=[inv1], 
                                    second_inversion=[inv2],
                                    t1map=qr1['t1'], r2smap=qr2map,
                                    uni=uni, 
                                    b1map=b1map['transformed_source'],
                                    b1scaling=1000.0,
                                    echo_times=[0.00199],
                                    inversion_times=[0.840, 2.370],
                                    flip_angles=[5.0, 6.0], inversion_TR=4.3,
                                    excitation_TR=[0.0072, 0.0072],
                                    N_excitations=180,
                                    save_data=True,overwrite=False,
                                    output_dir=proc_dir)
    
    else:
        # it turns out that B1 correction changes quite heavily the CSF / GM boundaries
        # this is handled later on, but stats might be quite different with or without B1 correction
        qr1 = nighres.intensity.mp2rage_t1_from_uni(uniform_image=uni,
                                        inversion_times=[0.840, 2.370],
                                        flip_angles=[5.0, 6.0], inversion_TR=4.3,
                                        excitation_TR=[0.0072, 0.0072],
                                        N_excitations=180,
                                        correct_B1=False, B1_map=None, 
                                        save_data=True,overwrite=False,
                                        file_name=subject+'_'+session+'_nob1',
                                        output_dir=proc_dir)
        
        qpd = nighres.intensity.mp2rageme_pd_mapping(first_inversion=[inv1], 
                                    second_inversion=[inv2],
                                    t1map=qr1nob1['t1'], r2smap=qr2map,
                                    uni=uni, 
                                    echo_times=[0.00199],
                                    inversion_times=[0.840, 2.370],
                                    flip_angles=[5.0, 6.0], inversion_TR=4.3,
                                    excitation_TR=[0.0072, 0.0072],
                                    N_excitations=180,
                                    save_data=True,overwrite=False,
                                    file_name=subject+'_'+session+'_nob1',
                                    output_dir=proc_dir)
    
    
    print("\n4. Skull stripping\n")
    
    inv2n4_file = os.path.join(proc_dir, nighres.utils._fname_4saving(rootfile=inv2,
                          suffix='n4'))
    
    if (not os.path.isfile(inv2n4_file)):
        #command = 'N4BiasFieldCorrection -d 3 -i '+inv2+' -o '+inv2n4_file
        #print(command)
        #try:
        #    subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        #except subprocess.CalledProcessError as e:
        #    msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        #    raise subprocess.CalledProcessError(msg)
        img = ants.image_read(inv2)
        inv2_n4 = ants.n4_bias_field_correction(img)
        ants.image_write(inv2_n4, inv2n4_file)
    
        # header issues with N4??
        img = nighres.io.load_volume(inv2)
        data = nighres.io.load_volume(inv2n4_file).get_fdata()
        img = nibabel.Nifti1Image(data,img.affine,img.header)
        nighres.io.save_volume(inv2n4_file, img)
            
    skull = nighres.brain.intensity_based_skullstripping(main_image=inv2n4_file, extra_image=qr1['t1'],
                        noise_model='exponential', 
                        skip_zero_values=True,
                        iterate=False, dilate_mask=0, topology_lut_dir=None,
                        save_data=True, overwrite=False, output_dir=proc_dir)
    
    brainmask = nighres.io.load_volume(skull['brain_mask']).get_fdata()
    
    #if b1ratio is not None:
    r1strip_file = proc_dir+subject+'_'+session+'_acq-mp2rage_mod-r1hz_orient-std_brain.nii.gz'
    if (not os.path.isfile(r1strip_file) and not os.path.isfile(final_qr1_file)):
        print("Mask qR1")
        r1 = nighres.io.load_volume(qr1['r1'])
        r1strip = nibabel.Nifti1Image(numpy.minimum(3.0,brainmask*r1.get_fdata()), r1.affine, r1.header)
    #    r1strip = nibabel.as_closest_canonical(r1strip)
        nighres.io.save_volume(r1strip_file, r1strip)
    
    r2strip_file = proc_dir+subject+'_'+session+'_acq-gre_mod-r2hz_orient-std_brain.nii.gz'
    if (not os.path.isfile(r2strip_file) and not os.path.isfile(final_qr2s_file)):
        print("Mask qR2*")
        r2s = nighres.io.load_volume(qr2map)
        r2strip = nibabel.Nifti1Image(numpy.minimum(200.0,brainmask*r2s.get_fdata()), r2s.affine, r2s.header)
    #    r2strip = nibabel.as_closest_canonical(r2strip)
        nighres.io.save_volume(r2strip_file, r2strip)
    
    # call N4
    pdn4_file = os.path.join(proc_dir, nighres.utils._fname_4saving(rootfile=qpd['pd'], suffix='n4'))
    if (not os.path.isfile(pdn4_file)):
        #command = 'N4BiasFieldCorrection -d 3 -i '+qpd['pd']+'  -x '+skull['brain_mask']+' -o '+pdn4_file
        #print(command)
        #try:
        #    subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        #except subprocess.CalledProcessError as e:
        #    msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        #    raise subprocess.CalledProcessError(msg)
        img = ants.image_read(qpd['pd'])
        msk = ants.image_read(skull['brain_mask'])
        pd_n4 = ants.n4_bias_field_correction(img, mask=msk)
        ants.image_write(pd_n4, pdn4_file)
    
    pdstrip_file = proc_dir+subject+'_'+session+'_acq-mp2rage_mod-qpd_orient-std_brain.nii.gz'
    if (not os.path.isfile(pdstrip_file) and not os.path.isfile(final_qpd_file)):
        print("Mask qPD")
        pd = nighres.io.load_volume(pdn4_file)
        pddata = pd.get_fdata()
        pdmean = numpy.mean(pddata[pddata>0])
        pdstrip = nibabel.Nifti1Image(numpy.minimum(8.0,brainmask*pddata/pdmean), pd.affine, pd.header)
    #    pdstrip = nibabel.as_closest_canonical(pdstrip)
        nighres.io.save_volume(pdstrip_file, pdstrip)
    
    
    print("\n5. QSM\n")

#    gre_phs_e2 = nighres.registration.apply_coordinate_mappings(gre_phs_e2, mapping['mapping'],
#                    interpolation="linear", padding="zero",
#                    save_data=True, overwrite=False, output_dir=proc_dir)['result']
#    gre_phs_e3 = nighres.registration.apply_coordinate_mappings(gre_phs_e3, mapping['mapping'],
#                    interpolation="linear", padding="zero",
#                    save_data=True, overwrite=False, output_dir=proc_dir)['result']
#    gre_phs_e4 = nighres.registration.apply_coordinate_mappings(gre_phs_e4, mapping['mapping'],
#                    interpolation="linear", padding="zero",
#                    save_data=True, overwrite=False, output_dir=proc_dir)['result']

    basic_header = nighres.io.load_volume(uni)

    # use the three longest echoes, compute QSM for each, take the median    
    phase2 = gre_phs_e2
    phase3 = gre_phs_e3
    phase4 = gre_phs_e4
    phases = [phase2,phase3,phase4]
    tes = [0.00722,0.01444,0.02323]

    # bring brain mask into GRE space
    brain_mask = nighres.registration.apply_coordinate_mappings(skull['brain_mask'], mapping['inverse'],
                                            interpolation="nearest", padding="zero",
                                            save_data=True, overwrite=False, output_dir=proc_dir)['result']

    for idx,phs in enumerate(phases):
        # check if already computed
        qsm_file = os.path.join(proc_dir, nighres.utils._fname_4saving(rootfile=phs,
                          suffix='tgv-qsm'))
    
        if (not os.path.isfile(qsm_file)):
            # here use the brain mask from earlier, but copy the header from phase
            # (the N4 step may have mangled the mask header...)
            qsm_mask_file = os.path.join(proc_dir, nighres.utils._fname_4saving(rootfile=phs,
                          suffix='tgv-mask'))
            #phs_img = nighres.io.load_volume(phs)
            mask_img = nibabel.Nifti1Image(nighres.io.load_volume(brain_mask).get_fdata(), basic_header.affine, basic_header.header)
            nighres.io.save_volume(qsm_mask_file, mask_img)
            mask = qsm_mask_file
            
            # rescale from [0,4095] to [0, 2PI]
            phs_file = os.path.join(proc_dir, nighres.utils._fname_4saving(rootfile=phs,
                          suffix='rad'))
            img = nighres.io.load_volume(phs)
            img = nibabel.Nifti1Image(img.get_fdata()/4095.0*2.0*math.pi,basic_header.affine,basic_header.header)
            nighres.io.save_volume(phs_file, img)
            
            # run tgv_qsm on each echo separately at first
            recon = 'tgv_qsm -p '+phs_file
            #recon = recon + ' -m '+mask['mask']
            recon = recon + ' -m '+mask
            recon = recon + ' -i '+str(1000)
            recon = recon + ' -f '+str(7)
            recon = recon + ' -t '+str(tes[idx])
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
            results = sorted(glob(proc_dir+'*_reconQSM*.nii.gz'))
            for result in results:
                os.rename(result, qsm_file)
    
    # recombine the three echos
    med_file = os.path.join(proc_dir, nighres.utils._fname_4saving(rootfile=phases[0],
                          suffix='tgv-qsm_med-img'))
    
    rng_file = os.path.join(proc_dir, nighres.utils._fname_4saving(rootfile=phases[0],
                          suffix='tgv-qsm_rng-img'))
    
    qsmstrip_file = proc_dir+subject+'_'+session+'_acq-gre_mod-qsm_orient-std_brain.nii.gz'
     
    if (not os.path.isfile(qsmstrip_file) and not os.path.isfile(final_qsm_file)):
        qsm2 = os.path.join(proc_dir, nighres.utils._fname_4saving(rootfile=phase2,
                          suffix='tgv-qsm'))
        qsm3 = os.path.join(proc_dir, nighres.utils._fname_4saving(rootfile=phase3,
                          suffix='tgv-qsm'))
        qsm4 = os.path.join(proc_dir, nighres.utils._fname_4saving(rootfile=phase4,
                          suffix='tgv-qsm'))
        
        data2 = nighres.io.load_volume(qsm2).get_data()
        data3 = nighres.io.load_volume(qsm3).get_data()
        data4 = nighres.io.load_volume(qsm4).get_data()
        
        data = numpy.stack((data2,data3,data4),axis=3)
        rng = numpy.max(data, axis=3) - numpy.min(data, axis=3)
        med = numpy.percentile(data, 50, axis=3)
        #avg = numpy.average(data, axis=3)
        
        brainmask = nighres.io.load_volume(brain_mask).get_fdata()
        qsmmask = scipy.ndimage.binary_erosion(brainmask, iterations=5)
        
        img = nighres.io.load_volume(qsm2)
        img = nibabel.Nifti1Image(qsmmask*numpy.maximum(-0.3,numpy.minimum(0.5,med)),img.affine,img.header)
        img = nibabel.as_closest_canonical(img)
    #    nighres.io.save_volume(med_file, img)
        nighres.io.save_volume(med_file, img)
    
        img = nibabel.Nifti1Image(rng,img.affine,img.header)
        nighres.io.save_volume(rng_file, img)
    
        # transform the QSM map to T1 space after processing
        qsmstrip_file = nighres.registration.apply_coordinate_mappings(med_file, mapping['mapping'],
                        interpolation="nearest", padding="zero",
                        save_data=True, overwrite=False, output_dir=proc_dir)['result']

    
    # copy to output folder with standardized naming convention
    print("6. Data transfer to target directory\n")
        
    if (not os.path.exists(res_dir)):
        os.makedirs(res_dir)
        
    if (not os.path.exists(final_qr1_file)):
        os.rename(r1strip_file,final_qr1_file)
    
    if (not os.path.exists(final_qr2s_file)):
        os.rename(r2strip_file,final_qr2s_file)
    
    if (not os.path.exists(final_qsm_file)):
        os.rename(qsmstrip_file,final_qsm_file)
    
    if (not os.path.exists(final_qpd_file)):
        os.rename(pdstrip_file,final_qpd_file)
    
