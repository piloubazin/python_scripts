
import nighres
import numpy
import math
import nibabel
import os
import shutil

subject='Ahead_brain_122017'
blockface = 'Ahead_brain_122017_blockface-image.nii.gz'

silver = 'Ahead_brain_122017_Bielschowsky-interpolated.nii.gz'
thio = 'Ahead_brain_122017_thionin-interpolated.nii.gz'
parv = 'Ahead_brain_122017_parvalbumin-interpolated.nii.gz'
calb = 'Ahead_brain_122017_calbindin-interpolated.nii.gz'
calret = 'Ahead_brain_122017_calretinin-interpolated.nii.gz'

output_dir = 'subcortex/'

ahead = nighres.data.download_AHEAD_template(data_dir=output_dir)

# co-register stains to AHEAD template
map_ahead = nighres.registration.embedded_antsreg(source_image=blockface, 
                    target_image=ahead['qr1'],
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
                    output_dir=output_dir,
                    file_name=subject+'_bf2ahead')

map_silver = nighres.registration.apply_coordinate_mappings(silver, map_ahead['mapping'], interpolation="nearest", padding="closest",
                        save_data=True, overwrite=True, output_dir=output_dir, file_name=subject+'_silver')

map_thio = nighres.registration.apply_coordinate_mappings(thio, map_ahead['mapping'], interpolation="nearest", padding="closest",
                        save_data=True, overwrite=True, output_dir=output_dir, file_name=subject+'_thio')

map_parv = nighres.registration.apply_coordinate_mappings(parv, map_ahead['mapping'], interpolation="nearest", padding="closest",
                        save_data=True, overwrite=True, output_dir=output_dir, file_name=subject+'_parv')

map_calb = nighres.registration.apply_coordinate_mappings(calb, map_ahead['mapping'], interpolation="nearest", padding="closest",
                        save_data=True, overwrite=False, output_dir=output_dir, file_name=subject+'_calb')

map_calret = nighres.registration.apply_coordinate_mappings(calret, map_ahead['mapping'], interpolation="nearest", padding="closest",
                        save_data=True, overwrite=False, output_dir=output_dir, file_name=subject+'_calret')

# transfer intensity priors
template_contrasts = [ahead['qr1'],ahead['qr2s'],ahead['qsm']]

target_contrasts = [blockface,silver,thio,parv,calb,calret]

aligned_contrasts = [map_ahead['transformed_source'],
                     map_silver['result'],
                     map_thio['result'],
                     map_parv['result'],
                     map_calb['result'],
                     map_calret['result']]

structures = ['str_hem-l','str_hem-r','stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r','rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
              'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r','vent_hem-3','vent_hem-4','amg_hem-l','amg_hem-r',
              'ic_hem-l','ic_hem-r','vta_hem-l','vta_hem-r','fx_hem-lr','pag_hem-l','pag_hem-r','ppn_hem-l','ppn_hem-r',
              'cl_hem-l','cl_hem-r']

n_contrasts = len(template_contrasts)
n_targets = len(target_contrasts)
n_structures = len(structures)

atlas = nighres.data.download_MASSP_atlas()

## first adapt the intensities
nighres.segmentation.conditional_shape_map_intensities(n_structures, n_contrasts, n_targets,
                      contrast_images=template_contrasts, 
                      target_images=aligned_contrasts,
                      shape_atlas_probas=atlas['spatial_probas'], 
                      shape_atlas_labels=atlas['spatial_labels'],
                      intensity_atlas_hist=atlas['histograms'],
                      skeleton_atlas_probas=atlas['skeleton_probas'], 
                      skeleton_atlas_labels=atlas['skeleton_labels'],
                      save_data=True, overwrite=True,
                      output_dir=output_dir, file_name=subject+'_insitu_contrasts')
                            
massp = nighres.parcellation.massp(target_images=target_contrasts,
                            map_to_target=map_ahead['inverse'], 
                            structures=n_structures, 
                            shape_atlas_probas=atlas['spatial_probas'], 
                            shape_atlas_labels=atlas['spatial_labels'],
                            intensity_atlas_hist=output_dir+subject+'_insitu_contrasts_cspmax-chist.nii.gz',
                            skeleton_atlas_probas=atlas['skeleton_probas'], 
                            skeleton_atlas_labels=atlas['skeleton_labels'],
                            max_iterations=120, max_difference=0.1,
                            save_data=True, file_name=subject+'_subcortex',
                            output_dir=output_dir, overwrite=True)

nighres.surface.parcellation_to_meshes(massp['max_label'], connectivity="18/6", 
                                         spacing = 0.0, smoothing=1.0,
                                         save_data=True, overwrite=True,
                                         output_dir=output_dir, file_name=subject+'_subcortex')