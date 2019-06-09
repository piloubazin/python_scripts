#!/usr/bin/python

import nibabel as nb
import numpy as np
import sys
import os
import argparse
import nighres as nh

parser = argparse.ArgumentParser()
parser.add_argument('command')
parser.add_argument('-o', '--outdir', nargs=1, help='output directory (default: input directory)')
parser.add_argument('imgfiles', nargs='+', help='input image files')
print(parser.parse_args(sys.argv))

imgfiles = parser.parse_args(sys.argv).imgfiles
outdirs = parser.parse_args(sys.argv).outdir
if outdirs==None : outdir = None
else : outdir = outdirs[0]

if len(imgfiles)==2:
	#nh.registration.embedded_syn(source_image=imgfiles[0], target_image=imgfiles[1], coarse_iterations=20, medium_iterations=0, fine_iterations=0,
	#				run_affine_first=True, cost_function='CrossCorrelation', interpolation='NearestNeighbor',
    #                save_data=True, output_dir=outdir,
    #                file_name=None)
	#nh.registration.embedded_syn(source_image=imgfiles[0], target_image=imgfiles[1], coarse_iterations=20, medium_iterations=0, fine_iterations=0)
	nh.laminar.laminar_iterative_smoothing(profile_surface_image=imgfiles[0], intensity_image=imgfiles[1], fwhm_mm=5,
                     save_data=True, output_dir=outdir,
                     file_name=None)


if len(imgfiles)==4:
    nh.cortex.cruise_cortex_extraction(init_image=imgfiles[0], wm_image=imgfiles[1], gm_image=imgfiles[2], csf_image=imgfiles[3], vd_image=None,
							data_weight=0.9, regularization_weight=0.1,
							max_iterations=800, normalize_probabilities=True,
							correct_wm_pv=True, wm_dropoff_dist=1.0,
							topology='wcs', topology_lut_dir=None,
							save_data=True, output_dir=outdir,
							file_name=None)