import nighres
import numpy
import nibabel
import glob
import os
import subprocess
mrtrix='/home/pilou/Software/Miniconda/bin/'

dwi_dir = '/home/public/StoppingAge-DWI-Ahead/Data/Dwi/'
qmr_dir = '/home/public/StoppingAge-DWI-Ahead/Data/qMRI-estimates/'
vsl_dir = '/home/pilou/Projects/VascularDiffusion/vessels/'
out_dir = '/home/pilou/Projects/VascularDiffusion/processing/'

qmr_myelin_files = glob.glob(qmr_dir+'sub-091/*-myelin.nii.gz')

overwrite=False

for qmr_file in qmr_myelin_files:
    # get subject id
    sub_id = qmr_file.replace(qmr_dir,'')[0:7]
    print('subject: '+sub_id)
    # find corresponding B0 map
    b0_file = glob.glob(dwi_dir+sub_id+'/preprocessed/mean_b0.nii.gz')
    if (len(b0_file)==0): print('mean B0 not found!')
    else:
        b0_file = b0_file[0]
        reg_dir = dwi_dir+sub_id+'/qmri_est/'
        tck_dir = dwi_dir+sub_id+'/tckgen/'
        wgt_dir = dwi_dir+sub_id+'/SIFT2-tckgen/'
        fod_dir = dwi_dir+sub_id+'/dwi2fod/'
        pre_dir = dwi_dir+sub_id+'/preprocessed/'
    
        coreg = nighres.registration.embedded_antsreg(source_image=qmr_file,target_image=b0_file,
                                                      run_rigid=True,run_affine=False,run_syn=False,
                                                      save_data=True,output_dir=reg_dir)
        
        vsl_probas = glob.glob(vsl_dir+sub_id+'*proba.nii.gz')
        if (len(vsl_probas)==0): print("missing vessel files!")
        else:
            for vsl_proba in vsl_probas:
                
                overwrite=False

                #vsl_propag = nighres.intensity.intensity_propagation(vsl_proba, combine='mean', distance_mm=1.275,
                #                                                        target='lower', scaling=0.5, 
                #                                                        save_data=True, overwrite=overwrite, output_dir=out_dir)
                
                trans = nighres.registration.apply_coordinate_mappings(image=vsl_proba,mapping1=coreg['mapping'],
                                                                        save_data=True, overwrite=overwrite, output_dir=out_dir)
                
                # make binarised version for selection
                mask = nighres.io.load_volume(trans['result'])
                mask = nibabel.Nifti1Image(mask.get_fdata()>0.5, mask.affine, mask.header)
                mask_file = out_dir+vsl_proba.replace(vsl_dir,'').replace('_mvf-proba.nii.gz','_mvf-proba_def-th0p5.nii.gz')
                nighres.io.save_volume(mask_file, mask)
                
                vsl_tract_file = out_dir+vsl_proba.replace(vsl_dir,'').replace('_mvf-proba.nii.gz', '_vdt-tract.tck')
                vsl_weight_file = out_dir+vsl_proba.replace(vsl_dir,'').replace('_mvf-proba.nii.gz', '_vdt-weights.txt')
                
                # start with recomputing the FODs over the entire brain? shouldn't be necessary...
                
                overwrite=False

                
                if ((os.path.exists(vsl_tract_file)) & (overwrite == False)):
                    print('\n Generate tracts: done (set overwrite to True to recompute)')
                else:
                    print('\n Generate tracts')
     
                    #command = mrtrix+'tckedit -include '+mask_file+' -tck_weights_in '+wgt_dir+'sift.txt '\
                    #                 +'-tck_weights_out '+vsl_weight_file\
                    #                    +' -force '+tck_dir+'tracks_combined.tck '+vsl_tract_file
        
                    command = mrtrix+'tckgen '+fod_dir+'wmfod.mif '+vsl_tract_file+' -seed_image '+mask_file\
                                    +' -mask '+pre_dir+sub_id+'_ses-3_dwi_mask.mif.gz'\
                                    +' -nthreads 4 -algorithm iFOD2 -select 100000 -force -cutoff 0.01'
                    
                    print(command+'\n')
                    try:
                        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
                    except subprocess.CalledProcessError as e:
                        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
                        raise subprocess.CalledProcessError(msg)
                
                
                vsl_map_file = out_dir+vsl_proba.replace(vsl_dir,'').replace('_mvf-proba.nii.gz', '_vdt-map.nii.gz')
                
                if ((os.path.exists(vsl_map_file)) & (overwrite == False)):
                    print('\n Generate tract map: done (set overwrite to True to recompute)')
                else:
                    print('\n Generate tract map')
 
                    #command = mrtrix+'tckmap -contrast scalar_map -image '+trans['result']+' -force '\
                    #              +'-tck_weights_in '+vsl_weight_file+' '+vsl_tract_file+' '+vsl_map_file
    
                    command = mrtrix+'tckmap -contrast scalar_map -image '+trans['result']+' -force '\
                                  +'-stat_vox max -stat_tck mean '\
                                  +'-template '+trans['result']+' '+vsl_tract_file+' '+vsl_map_file
    
                    print(command+'\n')
                    try:
                        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
                    except subprocess.CalledProcessError as e:
                        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
                        raise subprocess.CalledProcessError(msg)

                vsl_wgt_file = out_dir+vsl_proba.replace(vsl_dir,'').replace('_mvf-proba.nii.gz', '_vdt-wgt.txt')
                
                if ((os.path.exists(vsl_wgt_file)) & (overwrite == False)):
                    print('\n Generate tract weights: done (set overwrite to True to recompute)')
                else:
                    print('\n Generate tract weights')
 
                    #command = mrtrix+'tckmap -contrast scalar_map -image '+trans['result']+' -force '\
                    #              +'-tck_weights_in '+vsl_weight_file+' '+vsl_tract_file+' '+vsl_map_file
    
                    command = mrtrix+'tcksample -stat_tck mean -force '+vsl_tract_file+' '+trans['result']+' '+vsl_wgt_file
    
                    print(command+'\n')
                    try:
                        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
                    except subprocess.CalledProcessError as e:
                        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
                        raise subprocess.CalledProcessError(msg)

                overwrite = True

                vsl_sel_file = out_dir+vsl_proba.replace(vsl_dir,'').replace('_mvf-proba.nii.gz', '_vdt-sel.tck')
                
                if ((os.path.exists(vsl_sel_file)) & (overwrite == False)):
                    print('\n Generate tract selection: done (set overwrite to True to recompute)')
                else:
                    print('\n Generate tract selection')
 
                    #command = mrtrix+'tckmap -contrast scalar_map -image '+trans['result']+' -force '\
                    #              +'-tck_weights_in '+vsl_weight_file+' '+vsl_tract_file+' '+vsl_map_file
    
                    command = mrtrix+'tckedit -force -tck_weights_in '+vsl_wgt_file+' -minweight 0.1 -minlength 20.0 '\
                                        +vsl_tract_file+' '+vsl_sel_file
    
                    print(command+'\n')
                    try:
                        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
                    except subprocess.CalledProcessError as e:
                        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
                        raise subprocess.CalledProcessError(msg)
