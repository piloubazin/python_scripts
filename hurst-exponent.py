# small utility to compute DFA (detrended fluctuation analysis) and Hurst exponents

import nibabel as nb
import numpy as np
import sys
import os
import argparse
from scipy.fftpack import fft, ifft

parser = argparse.ArgumentParser()
parser.add_argument('command')
parser.add_argument('-o', '--outdir', nargs=1, help='output directory (default: input directory)')
parser.add_argument('-lf', '--low', nargs='?', default=0.01, help='lowest window frequency (default: 0.01 Hz)')
parser.add_argument('-hf', '--high', nargs='?', default=0.1, help='highest window frequency (default: 0.1 Hz)')
parser.add_argument('-rt', '--repetitiontime', nargs=1, help='fMRI repetition time (TR)')
parser.add_argument('-nw', '--nwindows', args='?', default=10, help='number of log-spaced windows to use (default: 10)')
parser.add_argument('-deg', '--degree', args='?', default=0, help='degree of the polynomial detrending (default: 0)')
parser.add_argument('imgfiles', nargs='+', help='input image files')
print(parser.parse_args(sys.argv))

# load the parameters
imgfiles = parser.parse_args(sys.argv).imgfiles
outdirs = parser.parse_args(sys.argv).outdir
low = float(parser.parse_args(sys.argv).low)
high = float(parser.parse_args(sys.argv).high)
tr = float(parser.parse_args(sys.argv).repetitiontime)
nw = int(parser.parse_args(sys.argv).nwindows)
deg = int(parser.parse_args(sys.argv).degree)
freq = 1.0/tr

for imgfile in imgfiles:
	# i/o info for files
	print("file to process: "+imgfile)
	if outdirs==None:
		outdir = os.path.dirname(imgfile)
		print("input/output directory: "+outdir)
	else:
		outdir=outdirs[0]
		print("output directory: "+outdir)

	imgdir = os.path.dirname(imgfile)

	imgname = os.path.basename(imgfile)
	basename = imgname.split(".")
	outname = outdir+basename[0]+"_dfa."+'.'.join(basename[1:])

	# loading the data
	print("opening file: "+imgfile)
	img = nb.load(imgfile)
	data = img.get_data()
	affine = img.affine #or img.get_affine(), which I think is being deprecated?
	header = img.header

	length = np.shape(data)[3]

	print("defining the frequency window")
	nextpowerof2 = np.ceil(np.log2(length))
	padded = int(np.power(2, nextpowerof2))
	
	#print("freq:",freq,"lf:",low,"hf:",high,"length:",padded)
	if (low >= freq/2):
		lowid = int(padded/2)
	else:
		lowid = int(np.ceil(low*padded*tr))
	
	if (high >= freq/2):
		highid = int(padded/2)
	else:
		highid = int(np.floor(high*padded*tr))
	
	frequencymask = np.zeros(padded)
	frequencymask[lowid+1:highid+1] = 1
	frequencymask[padded-highid:padded-lowid] = 1
	#print(lowid,highid)
	#print(frequencymask)
	
	print("removing the mean")
	datamean = data.mean(3,keepdims=True)
	datamean = np.repeat(datamean,length,axis=3)
	
	data = np.pad(data-datamean,((0,0),(0,0),(0,0),(0,padded-length)),'constant') 
	#print("shape: ",np.shape(data))
	
	print("filtering")
	data = np.fft.fft(data,axis=3)
	#print("shape: ",np.shape(data))
	data[:,:,:,np.nonzero(frequencymask==0)] = 0
	data = np.real(np.fft.ifft(data,axis=3))
	
	data = data[:,:,:,0:length]+datamean
		
	print("saving to: "+outname)
	outimg = nb.Nifti1Image(data, affine, header)
	outimg.header['cal_min'] = np.min(data)
	outimg.header['cal_max'] = np.max(data)
	outimg.to_filename(outname)
