# small python script to handle the various dwi processing steps.

import nighres
import numpy
import nibabel
import subprocess
import os
import glob
from scripts.total_readout_time_dcm import total_readout_time_dcm
import shutil

# subject ID
subject = 'sub-xxx'

# where is the DICOM MR input folder located
dicom_dir = '/home/pilou/Projects/Diffusion/ahead-DWI-piloting-Luka-Matthan/feb19/20190211/diffusion_pilot3T_20190211_diffusion_pilot3T_20190211/20190211_1.3.46.670589.11.78004.5.0.31716.2019021112061994000/MR'
dicom_img = '00701_dwi64_1stHALF_1.28mm_b1600/09160.dcm'
dicom_tp = '00601_dwi_1.28mm_b0_TOPUP/00972.dcm'

# where we dump all intermediate files
work_dir = '/home/pilou/Projects/Diffusion/ahead-DWI-piloting-Luka-Matthan/feb19/20190211/preprocessing'

# where we store the result
output_dir = '/home/pilou/Projects/Diffusion/ahead-DWI-piloting-Luka-Matthan/feb19/20190211/dwi'


if not os.path.exists(work_dir):
    print('Create temporary work directory')
    os.mkdir(work_dir)

## 1. Convert dat from DICOM to Nifti using dcm2niix
convert_dir = '/dcm2niix'
if (os.path.exists(work_dir+convert_dir)):
    print('Conversion: done (delete dir to recompute)')
else:
    os.mkdir(work_dir+convert_dir)
    
    print('\n1. Convert DICOM into Nifti files')
    command = 'dcm2niix -z y -f S%s_%d -o '+work_dir+convert_dir+' '+dicom_dir
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)

    ## 2. Convert to MRtrix3 format: only files with associated bval,bvec
    print('\n2. Convert Nifti files into Mrtrix3 files')
    bval_files = glob.glob(work_dir+convert_dir+'/*.bval')
    for bval_file in bval_files:
        bvec_file = bval_file.replace('.bval','.bvec')
        data_file = bval_file.replace('.bval','.nii.gz')
        mif_file = bval_file.replace('.bval','.mif.gz')
        
        command = 'mrconvert -fslgrad '+bvec_file+' '+bval_file+' '+data_file+' '+mif_file
    
        print(command)
        try:
            subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
            raise subprocess.CalledProcessError(msg)

    print('\n CHECK: are there any files you need to remove from among these: '+str(glob.glob(work_dir+convert_dir+'/*.mif.gz')))

    proceed = input('proceed?(y/n): ')
    if proceed is not 'y':
        sys.exit(1)
    
## before continuing: manually remove .mif.gz files if some of the diffusion files are not used

## 3. Concatenate all dwi images into a single one
concat_dir = '/concat'
if (os.path.exists(work_dir+concat_dir)):
    print('Data combination: done (delete dir to recompute)')
else:
    print('\n3. Concatenate dwi and topup files into a single file each')
    os.mkdir(work_dir+concat_dir)
    mif_files = glob.glob(work_dir+convert_dir+'/*.mif.gz')
    
    dwi_files = []
    topup_files = []
    for mif_file in mif_files:
        if mif_file.find('TOPUP')>-1 or mif_file.find('topup')>-1:
            topup_files.append(mif_file)
        else:
            dwi_files.append(mif_file)
    
    print('\nmerge '+str(len(dwi_files))+' dwi files')
    command = 'mrcat -axis 3 '
    for dwi_file in dwi_files:
        command = command+dwi_file+' '
    command = command+work_dir+concat_dir+'/dti_raw.mif.gz'
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)

    print('\nmerge '+str(len(topup_files))+' topup files')
    command = 'mrcat -axis 3 '
    for topup_file in topup_files:
        command = command+topup_file+' '
    command = command+work_dir+concat_dir+'/dti_topup_raw.mif.gz'
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)

## 4. Switching to Luka Liebrand's customized script (recoded in python for general use)
denoise_dir = '/denoise'
if (os.path.exists(work_dir+denoise_dir)):
    print('Denoising: done (delete dir to recompute)')
else:
    print('\n4. Denoise raw DWI files')
    os.mkdir(work_dir+denoise_dir)

    command = 'dwidenoise -nthreads 4 -extent 3,3,3'
    command = command +' -noise '+work_dir+denoise_dir+'/dti_noisemap.mif.gz'
    command = command +' '+work_dir+concat_dir+'/dti_raw.mif.gz'
    command = command +' '+work_dir+denoise_dir+'/dti_denoised.mif.gz'

    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)

    command = 'mrcalc '+work_dir+concat_dir+'/dti_raw.mif.gz'
    command = command +' '+work_dir+denoise_dir+'/dti_denoised.mif.gz'
    command = command +' -subtract '+work_dir+denoise_dir+'/dti_residuals.mif.gz'
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)


denoise_tp_dir = '/denoise_tp'
if (os.path.exists(work_dir+denoise_tp_dir)):
    print('Denoising topup: done (delete dir to recompute)')
else:
    print('\n4. Denoise raw topup DWI files')
    os.mkdir(work_dir+denoise_tp_dir)

    command = 'dwidenoise -nthreads 4 -extent 3,3,3'
    command = command +' -noise '+work_dir+denoise_tp_dir+'/dti_tp_noisemap.mif.gz'
    command = command +' '+work_dir+concat_dir+'/dti_topup_raw.mif.gz'
    command = command +' '+work_dir+denoise_tp_dir+'/dti_tp_denoised.mif.gz'

    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)

    command = 'mrcalc '+work_dir+concat_dir+'/dti_topup_raw.mif.gz'
    command = command +' '+work_dir+denoise_tp_dir+'/dti_tp_denoised.mif.gz'
    command = command +' -subtract '+work_dir+denoise_tp_dir+'/dti_tp_residuals.mif.gz'
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)


degibbs_dir = '/degibbs'
if (os.path.exists(work_dir+degibbs_dir)):
    print('De-Gibbs: done (delete dir to recompute)')
else:
    print('\n5a. Removing Gibbs ringing artefact from denoised DWI files')
    os.mkdir(work_dir+degibbs_dir)

    ## here we adapt Matthan Caan's de-Gibbs script ##
    # 1. Convert back to Nifti
    command = 'mrconvert -export_grad_fsl '+work_dir+degibbs_dir+'/dti_denoised.bvec'
    command = command +' '+work_dir+degibbs_dir+'/dti_denoised.bval'
    command = command +' '+work_dir+denoise_dir+'/dti_denoised.mif.gz'
    command = command +' '+work_dir+degibbs_dir+'/dti_denoised.nii.gz'
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)
    
    # 2. Do the processing
    print('Compute second derivatives for B0 images')
    infile = open(work_dir+degibbs_dir+'/dti_denoised.bval', 'r')
    bvals = [float(w) for w in infile.read().split()]
    infile.close()
    
    img = nighres.io.load_volume(work_dir+degibbs_dir+'/dti_denoised.nii.gz')
    data = img.get_data()

    # only for B0 values (e.g.bval < 100)
    for idx,bval in enumerate(bvals):
        if bval<100:
            print('found B0: '+str(bval))
            bzero = data[:,:,:,idx]
    
            deriv = numpy.gradient(numpy.gradient(bzero,axis=0),axis=0)+numpy.gradient(numpy.gradient(bzero,axis=1),axis=1)
            deriv[deriv<0] = 0
    
            bzero = bzero + deriv
            data[:,:,:,idx] = bzero
    
    img = nibabel.Nifti1Image(data, affine=img.affine,header=img.header)
    nighres.io.save_volume(work_dir+degibbs_dir+'/dti_unring.nii.gz', img)
    
    # 3. Convert back not needed: next step requires .nii.gz


# alternatively: use the MRtrix3 tools?
degibbs2_dir = '/degibbs2'
if (os.path.exists(work_dir+degibbs2_dir)):
    print('De-Gibbs ringing by MRtrix: done (delete dir to recompute)')
else:
    print('\n5b. De-Gibbs denoised DWI files (MRtrix version)')
    os.mkdir(work_dir+degibbs2_dir)

    command = 'mrdegibbs -nthreads 4 -axes 0,1'
    command = command +' '+work_dir+denoise_dir+'/dti_denoised.mif.gz'
    command = command +' '+work_dir+degibbs2_dir+'/dti_unring2.mif.gz'

    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)

    command = 'mrconvert -export_grad_fsl '+work_dir+degibbs2_dir+'/dti_unring2.bvec'
    command = command +' '+work_dir+degibbs2_dir+'/dti_unring2.bval'
    command = command +' '+work_dir+degibbs2_dir+'/dti_unring2.mif.gz'
    command = command +' '+work_dir+degibbs2_dir+'/dti_unring2.nii.gz'
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)


degibbs_tp_dir = '/degibbs_tp'
if (os.path.exists(work_dir+degibbs_tp_dir)):
    print('De-Gibbs topup: done (delete dir to recompute)')
else:
    print('\n5c. Removing Gibbs ringing artefact from denoised topup DWI files')
    os.mkdir(work_dir+degibbs_tp_dir)

    ## here we adapt Matthan Caan's de-Gibbs script ##
    # 1. Convert back to Nifti
    command = 'mrconvert -export_grad_fsl '+work_dir+degibbs_tp_dir+'/dti_tp_denoised.bvec'
    command = command +' '+work_dir+degibbs_tp_dir+'/dti_tp_denoised.bval'
    command = command +' '+work_dir+denoise_tp_dir+'/dti_tp_denoised.mif.gz'
    command = command +' '+work_dir+degibbs_tp_dir+'/dti_tp_denoised.nii.gz'
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)
    
    # 2. Do the processing
    print('Compute second derivatives for B0 images')
    infile = open(work_dir+degibbs_tp_dir+'/dti_tp_denoised.bval', 'r')
    bvals = [float(w) for w in infile.read().split()]
    infile.close()
    
    img = nighres.io.load_volume(work_dir+degibbs_tp_dir+'/dti_tp_denoised.nii.gz')
    data = img.get_data()

    # only for B0 values (e.g.bval < 100)
    for idx,bval in enumerate(bvals):
        if bval<100:
            print('found B0: '+str(bval))
            bzero = data[:,:,:,idx]
    
            deriv = numpy.gradient(numpy.gradient(bzero,axis=0),axis=0)+numpy.gradient(numpy.gradient(bzero,axis=1),axis=1)
            deriv[deriv<0] = 0
    
            bzero = bzero + deriv
            data[:,:,:,idx] = bzero
    
    img = nibabel.Nifti1Image(data, affine=img.affine,header=img.header)
    nighres.io.save_volume(work_dir+degibbs_tp_dir+'/dti_tp_unring.nii.gz', img)
    
    # 3. Convert back not needed: next step requires .nii.gz

# alternatively: use the MRtrix3 tools?
degibbs2_tp_dir = '/degibbs2_tp'
if (os.path.exists(work_dir+degibbs2_tp_dir)):
    print('De-Gibbs topup ringing by MRtrix: done (delete dir to recompute)')
else:
    print('\n5d. De-Gibbs denoised topup DWI files (MRtrix version)')
    os.mkdir(work_dir+degibbs2_tp_dir)

    command = 'mrdegibbs -nthreads 4 -axes 0,1'
    command = command +' '+work_dir+denoise_tp_dir+'/dti_tp_denoised.mif.gz'
    command = command +' '+work_dir+degibbs2_tp_dir+'/dti_tp_unring2.mif.gz'

    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)

    command = 'mrconvert -export_grad_fsl '+work_dir+degibbs2_tp_dir+'/dti_tp_unring2.bvec'
    command = command +' '+work_dir+degibbs2_tp_dir+'/dti_tp_unring2.bval'
    command = command +' '+work_dir+degibbs2_tp_dir+'/dti_tp_unring2.mif.gz'
    command = command +' '+work_dir+degibbs2_tp_dir+'/dti_tp_unring2.nii.gz'
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)


topup_dir = '/topup'
if (os.path.exists(work_dir+topup_dir)):
    print('Topup: done (delete dir to recompute)')
else:
    print('\n6. Top-up geometric distortion correction')
    os.mkdir(work_dir+topup_dir)

    # top-up requires an even number of slices: pad if necessary
    img = nighres.io.load_volume(work_dir+degibbs_dir+'/dti_unring.nii.gz')
    data = img.get_data()

    if (data.shape[2]%2==1):
        padded = numpy.zeros((data.shape[0],data.shape[1],data.shape[2]+1,data.shape[3]))
        padded[0:data.shape[0],0:data.shape[1],0:data.shape[2],0:data.shape[3]] = data
    else:
        padded = data
        
    img = nibabel.Nifti1Image(padded, affine=img.affine,header=img.header)
    nighres.io.save_volume(work_dir+topup_dir+'/dti.nii.gz', img)
    
    shutil.copy(work_dir+degibbs_dir+'/dti_denoised.bval',work_dir+topup_dir+'/dti.bval')
    shutil.copy(work_dir+degibbs_dir+'/dti_denoised.bvec',work_dir+topup_dir+'/dti.bvec')
    
    img = nighres.io.load_volume(work_dir+degibbs_tp_dir+'/dti_tp_unring.nii.gz')
    data = img.get_data()

    if (data.shape[2]%2==1):
        padded = numpy.zeros((data.shape[0],data.shape[1],data.shape[2]+1,data.shape[3]))
        padded[0:data.shape[0],0:data.shape[1],0:data.shape[2],0:data.shape[3]] = data
    else:
        padded = data
        
    img = nibabel.Nifti1Image(padded, affine=img.affine,header=img.header)
    nighres.io.save_volume(work_dir+topup_dir+'/dti_tp.nii.gz', img)
    
    shutil.copy(work_dir+degibbs_tp_dir+'/dti_tp_denoised.bval',work_dir+topup_dir+'/dti_tp.bval')
    shutil.copy(work_dir+degibbs_tp_dir+'/dti_tp_denoised.bvec',work_dir+topup_dir+'/dti_tp.bvec')
    
    # run the fsl script
    # 1. (extract total readout time from DICOM files)
    totreadout_main = total_readout_time_dcm(dicom_dir+'/'+dicom_img)
    totreadout_topup = total_readout_time_dcm(dicom_dir+'/'+dicom_tp)
    
    # 2. (merge all B0 from normal and topup acquisitions)
    command = 'dwiextract -bzero'
    command = command +' -fslgrad '+work_dir+topup_dir+'/dti.bvec '+work_dir+topup_dir+'/dti.bval'
    command = command +' '+work_dir+topup_dir+'/dti.nii.gz'
    command = command +' '+work_dir+topup_dir+'/b0.nii.gz'
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)
    
    command = 'dwiextract -bzero'
    command = command +' -fslgrad '+work_dir+topup_dir+'/dti_tp.bvec '+work_dir+topup_dir+'/dti_tp.bval'
    command = command +' '+work_dir+topup_dir+'/dti_tp.nii.gz'
    command = command +' '+work_dir+topup_dir+'/b0_tp.nii.gz'
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)
    
    command = 'mrcat -axis 3 '
    command = command +' '+work_dir+topup_dir+'/b0.nii.gz'
    command = command +' '+work_dir+topup_dir+'/b0_tp.nii.gz'
    command = command +' '+work_dir+topup_dir+'/b0_all.nii.gz'
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)
    
    # 3. (build the corresponding acquisition parameter file)
    acqfile = open(work_dir+topup_dir+'/acqparams.txt', 'w')
    b0_img = nighres.io.load_volume(work_dir+topup_dir+'/b0.nii.gz')
    b0_tp_img = nighres.io.load_volume(work_dir+topup_dir+'/b0_tp.nii.gz')
    for i in range(b0_img.header.get_data_shape()[3]):
        acqfile.write("0 1 0 "+str(totreadout_main)+"\n")
    for i in range(b0_tp_img.header.get_data_shape()[3]):
        acqfile.write("0 -1 0 "+str(totreadout_topup)+"\n")
    acqfile.close()
    
    idxfile = open(work_dir+topup_dir+'/index.txt', 'w')
    dti_img = nighres.io.load_volume(work_dir+topup_dir+'/dti.nii.gz')
    for i in range(dti_img.header.get_data_shape()[3]):
        idxfile.write("1  ")
    idxfile.close()
        
    # 4. (run topup)
    command = 'topup --imain='+work_dir+topup_dir+'/b0_all.nii.gz'
    command = command +' --datain='+work_dir+topup_dir+'/acqparams.txt'
    command = command +' --config=scripts/b02b0_hires.cnf'
    command = command +' --out='+work_dir+topup_dir+'/topup'
    command = command +' --fout='+work_dir+topup_dir+'/tfield'
    command = command +' --iout='+work_dir+topup_dir+'/topuped'
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)
    
    # 5. (average A and P DWIs)
    command = 'fslmaths '+work_dir+topup_dir+'/topuped -Tmean '+work_dir+topup_dir+'/meanB0_AP'
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)
    
    # 6. (mask with BET)
    command = 'bet '+work_dir+topup_dir+'/meanB0_AP '+work_dir+topup_dir+'/bet_meanB0_AP -f 0.1 -m'
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)
    


eddy_dir = '/eddy'
if (os.path.exists(work_dir+eddy_dir)):
    print('Eddy: done (delete dir to recompute)')
else:
    print('\n7. Eddy current distortion correction')
    os.mkdir(work_dir+eddy_dir)

    # 7. (eddy)
    command = 'eddy_openmp --imain='+work_dir+topup_dir+'/dti.nii.gz'
    command = command +' --mask='+work_dir+topup_dir+'/bet_meanB0_AP_mask.nii.gz'
    command = command +' --acqp='+work_dir+topup_dir+'/acqparams.txt'
    command = command +' --index='+work_dir+topup_dir+'/index.txt'
    command = command +' --bvecs='+work_dir+topup_dir+'/dti.bvec'
    command = command +' --bvals='+work_dir+topup_dir+'/dti.bval'
    command = command +' --topup='+work_dir+topup_dir+'/topup'
    command = command +' --repol'
    command = command +' --out='+work_dir+eddy_dir+'/eddy_corr_data'
    command = command +' --nvoxhp=4000 --ol_nstd=3 --slm=linear --mb=2 --fwhm=10,0,0,0,0'

    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)
    
    
if (os.path.exists(output_dir) and os.path.exists(output_dir+'/'+subject+'_dwi-preproc.mif.gz')):
    print('Final DWI image: done (delete to recompute)')
else:
    print('\n8. Final image creation')
    os.mkdir(output_dir)

    command = 'mrconvert -fslgrad '+work_dir+topup_dir+'/dti.bvec'
    command = command +' '+work_dir+topup_dir+'/dti.bval'
    command = command +' '+work_dir+eddy_dir+'/eddy_corr_data.nii.gz'
    command = command +' '+output_dir+'/'+subject+'_dwi-preproc.mif.gz'
    
    print(command)
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
        raise subprocess.CalledProcessError(msg)

    