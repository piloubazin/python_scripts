#!/usr/bin/python

# adapted from SBL bandpass filtering toolobox in Matlab

import nibabel as nb
import numpy as np
import sys
import os
import argparse
from scipy.fftpack import fft, ifft

parser = argparse.ArgumentParser()
parser.add_argument('command')
parser.add_argument('-o', '--outdir', nargs=1, help='output directory (default: input directory)')
parser.add_argument('-f', '--fmri', nargs='*', metavar='sub#:run1,run2...', help='list of fMRI time series (4D) per run per subject')
parser.add_argument('-m', '--mask', nargs='*', metavar='sub#:run1,run2...', help='list of masks (3D) per run per subject')
parser.add_argument('-t', '--transform', nargs='*', metavar='sub#:run1,run2...', help='list of transformations (coordinate mappings) to group space per run per subject')
parser.add_argument('-mc', '--meancutoff', nargs='?', default=0.0001, help='mask out voxels with (normalized) mean over time below cutoff')
parser.add_argument('-gc', '--groupcutoff', nargs='?', default=0.5, help='discard voxels with fewer subjects than this fraction of the group')
parser.add_argument('-sf', '--skip_first', nargs='?', default=0, help='skip time frames at the beginning of each run')
parser.add_argument('-sl', '--skip_last', nargs='?', default=0, help='skip time frames at the end of each run')
parser.add_argument('-nr', '--nruns', nargs='?', default=1, help='number of runs per subject')
print(parser.parse_args(sys.argv))

# load the parameters
fmrifiles = parser.parse_args(sys.argv).fmri
maskfiles = parser.parse_args(sys.argv).mask
transformfiles = parser.parse_args(sys.argv).transform
outdir = parser.parse_args(sys.argv).outdir
if outdir==None:
    outdir = os.getcwd()
    print("input/output directory: "+outdir)
else:
    outdir = outdir[0]
    print("output directory: "+outdir)
if (outdir[-1]!='/') : outdir +='/'

meancutoff = float(parser.parse_args(sys.argv).meancutoff)
groupcutoff = float(parser.parse_args(sys.argv).groupcutoff)
skip_first = int(parser.parse_args(sys.argv).skip_first)
skip_last = int(parser.parse_args(sys.argv).skip_last)
nruns = int(parser.parse_args(sys.argv).nruns)
nsubjects = int(len(fmrifiles)/nruns)


basename = os.path.basename(fmrifiles[0])
basename = basename.split(".")

print('ISC: '+str(nsubjects)+" subjects x "+str(nruns)+" runs, (mc,gc,sf,sl): ",meancutoff,groupcutoff,skip_first,skip_last)

# for debugging / information
debug = True
verbose = True

# for saving intermediate data
savezscores = True
saveavgfmri = True

# convenience labels
X=0
Y=1
Z=2
T=3

def main() :

    # 1. build a global average of all subjects (removing voxels with low intensity, missing values, and z-scoring the rest)
    (avgfmri, avgcount, avgmask, p0, pM, runtimes) = build_average(fmrifiles, maskfiles, transformfiles, 
                                                                    nsubjects, nruns, skip_first, skip_last, 
                                                                    meancutoff, groupcutoff)
    
    # 2. for each subject, remove data from the average and compute the correlation
    correlation = compute_subject_correlations(fmrifiles, maskfiles, transformfiles, 
                                               nsubjects, nruns, skip_first, skip_last, 
                                               avgfmri, avgcount, avgmask, p0, pM, runtimes)
    
    # export results (resampled to avg space)
    transform = nb.load(transformfiles[0])
    nx = transform.header.get_data_shape()[X]
    ny = transform.header.get_data_shape()[Y]
    nz = transform.header.get_data_shape()[Z]
    resampled = np.zeros((nx,ny,nz,nsubjects))
    for x in xrange(p0[X],pM[X]):
        for y in xrange(p0[Y],pM[Y]):
            for z in xrange(p0[Z],pM[Z]):
                if (avgmask[x,y,z]>0) :
                    for s in xrange(nsubjects):
                        resampled[x,y,z,s] = correlation[x-p0[X],y-p0[Y],z-p0[Z],s]
                        
    outname = outdir+basename[0]+"_isc"+str(nsubjects)+"x"+str(nruns)+"_corr."+".".join(basename[1:])
    print("saving to: "+outname)
    outimg = nb.Nifti1Image(resampled, transform.affine, transform.header)
    outimg.header['cal_min'] = np.min(resampled)
    outimg.header['cal_max'] = np.max(resampled)
    outimg.to_filename(outname)
    
    # 3. additional statistics?

    return


def compute_subject_correlations(fmrifiles, maskfiles, transformfiles, nsubjects, nruns, skip_first, skip_last, avgfmri, avgcount, avgmask, p0, pM, runtimes) :

    # test for consistency
    # here we assume the same numbers of subjects x runs for all data
    # (different versions could be made for simpler cases)
    if not ( len(fmrifiles) == len(maskfiles) and len(fmrifiles) == len(transformfiles) and len(fmrifiles) == nsubjects*nruns):
        print("!different number of masks, transform and fmri data than specified subjects and runs!")
        return

    # for each subject, remove from average and correlate
    correlation = np.zeros((pM[X]-p0[X],pM[Y]-p0[Y],pM[Z]-p0[Z],nsubjects))
    for s in xrange(nsubjects) :
        # keep sums in case of multiple runs per subject
        
        covar = np.zeros((pM[X]-p0[X],pM[Y]-p0[Y],pM[Z]-p0[Z]))
        var_sub = np.zeros((pM[X]-p0[X],pM[Y]-p0[Y],pM[Z]-p0[Z]))
        var_avg = np.zeros((pM[X]-p0[X],pM[Y]-p0[Y],pM[Z]-p0[Z]))
        samples = np.zeros((pM[X]-p0[X],pM[Y]-p0[Y],pM[Z]-p0[Z]))
        for r in xrange(nruns) :
            subjectrun = nb.load(fmrifiles[r+s*nruns])
            fmri = subjectrun.get_data()
            transform = nb.load(transformfiles[r+s*nruns]).get_data()
            mask = nb.load(maskfiles[r+s*nruns]).get_data()
        
            fmri = build_zscore(fmri, mask, transform, skip_first, skip_last, meancutoff)
            
            # remove from global average
            for x in xrange(p0[X],pM[X]):
                for y in xrange(p0[Y],pM[Y]):
                    for z in xrange(p0[Z],pM[Z]):
                        if (avgmask[x,y,z]>0) :
                            xp = int(np.rint(transform[x,y,z,X]))
                            yp = int(np.rint(transform[x,y,z,Y]))
                            zp = int(np.rint(transform[x,y,z,Z]))
                        
                            if (mask[xp,yp,zp]>0) : # only where you have valid data
                                for t in xrange(skip_first, fmri.shape[T]-skip_last) :
                                    # remove normalization
                                    avgdata = avgfmri[x-p0[X],y-p0[Y],z-p0[Z],runtimes[r]+t-skip_first]*avgcount[x-p0[X],y-p0[Y],z-p0[Z],r]
                                    # remove data value
                                    avgdata -= fmri[xp,yp,zp,t]
                                    # re-normalize
                                    avgdata /= (avgcount[x-p0[X],y-p0[Y],z-p0[Z],r]-1)
                                    # update variances and covariance
                                    covar[x-p0[X],y-p0[Y],z-p0[Z]] += fmri[xp,yp,zp,t]*avgdata
                                    var_sub[x-p0[X],y-p0[Y],z-p0[Z]] += fmri[xp,yp,zp,t]*fmri[xp,yp,zp,t]
                                    var_avg[x-p0[X],y-p0[Y],z-p0[Z]] += avgdata*avgdata
                                    samples[x-p0[X],y-p0[Y],z-p0[Z]] += 1

        # compute the final correlations
        for x in xrange(p0[X],pM[X]):
            for y in xrange(p0[Y],pM[Y]):
                for z in xrange(p0[Z],pM[Z]):
                    if (avgmask[x,y,z]>0) :
                        correlation[x-p0[X],y-p0[Y],z-p0[Z],s] = covar[x-p0[X],y-p0[Y],z-p0[Z]] \
                                                                /np.sqrt(var_sub[x-p0[X],y-p0[Y],z-p0[Z]]) \
                                                                /np.sqrt(var_avg[x-p0[X],y-p0[Y],z-p0[Z]])
                                               
    return correlation                       

def build_average( fmrifiles, maskfiles, transformfiles, nsubjects, nruns, skip_first, skip_last, meancutoff, groupcutoff) :
    
    # test for consistency
    # here we assume the same numbers of subjects x runs for all data
    # (different versions could be made for simpler cases)
    if not ( len(fmrifiles) == len(maskfiles) and len(fmrifiles) == len(transformfiles) and len(fmrifiles) == nsubjects*nruns):
        print("!different number of masks, transform and fmri data than specified subjects and runs!")
        return
        
    # init from first subject
    if debug : print("opening file: "+maskfiles[0])
    mask = nb.load(maskfiles[0]).get_data()
    
    if debug : print("opening file: "+transformfiles[0])
    transform = nb.load(transformfiles[0]).get_data()
    
    # define common space: every non-zero element of the mask is used
    nx = transform.shape[X]
    ny = transform.shape[Y]
    nz = transform.shape[Z]
    
    avgmask = np.zeros(transform.shape[X:T])
    for n in xrange(0,nsubjects*nruns) :
        if debug : print("subject ",n)
        transform = nb.load(transformfiles[n]).get_data()
        mask = nb.load(maskfiles[n]).get_data()

        for x in xrange(nx):
            for y in xrange(ny):
                for z in xrange(nz):
                    xp = int(np.rint(transform[x,y,z,X]))
                    yp = int(np.rint(transform[x,y,z,Y]))
                    zp = int(np.rint(transform[x,y,z,Z]))
                    avgmask[x,y,z] += mask[xp,yp,zp]
                    
    # find mask boundaries
    indices = np.nonzero(avgmask)
    p0 = (int(np.min(indices[X])),int(np.min(indices[Y])),int(np.min(indices[Z])))
    pM = (int(np.max(indices[X])+1),int(np.max(indices[Y])+1),int(np.max(indices[Z])+1))
    print("analysis bounding box (global space): ", p0, pM)
    
    # create average fmri only within the mask, concatenate the runs
    runtimes = range(nruns+1)
    # first run time is zero, simpler for later computations
    runtimes[0] = 0
    for n in xrange(nruns) :
        fmri_time = nb.load(fmrifiles[n]).header.get_data_shape()[T]
        runtimes[n+1] = int(runtimes[n] + fmri_time-skip_first-skip_last)
        
    avgfmri = np.zeros((pM[X]-p0[X],pM[Y]-p0[Y],pM[Z]-p0[Z],runtimes[nruns]))
    # to count the number of subjects per run that contribute to each voxel
    avgcount = np.zeros((pM[X]-p0[X],pM[Y]-p0[Y],pM[Z]-p0[Z],nruns))
        
    # pull data from each time series
    # check for % missing data (below threshold) and mask out the bad ones
    # use the rest to build Z-scores
    # note: for simplicity, the number of missing data points is the maximum over runs per voxel
    for s in xrange(nsubjects) :
        for r in xrange(nruns) :
            subjectrun = nb.load(fmrifiles[r+s*nruns])
            fmri = subjectrun.get_data()
            transform = nb.load(transformfiles[r+s*nruns]).get_data()
            mask = nb.load(maskfiles[r+s*nruns]).get_data()
        
            fmri = build_zscore(fmri, mask, transform, skip_first, skip_last, meancutoff)
            # option: save the masked, z-scored time series?
            if savezscores :
                imgname = os.path.basename(fmrifiles[r+s*nruns])
                basename = imgname.split(".")
                outname = outdir+basename[0]+"_zscored."+'.'.join(basename[1:])
                print("saving to: "+outname)
                outimg = nb.Nifti1Image(fmri, subjectrun.affine, subjectrun.header)
                outimg.header['cal_min'] = np.min(fmri)
                outimg.header['cal_max'] = np.max(fmri)
                outimg.to_filename(outname)	
    
            # build the global average
            for x in xrange(p0[X],pM[X]):
                for y in xrange(p0[Y],pM[Y]):
                    for z in xrange(p0[Z],pM[Z]):
                        xp = int(np.rint(transform[x,y,z,X]))
                        yp = int(np.rint(transform[x,y,z,Y]))
                        zp = int(np.rint(transform[x,y,z,Z]))
                        
                        if (mask[xp,yp,zp]>0) :
                            for t in xrange(skip_first, fmri.shape[T]-skip_last) :
                                avgfmri[x-p0[X],y-p0[Y],z-p0[Z],runtimes[r]+t-skip_first] += fmri[xp,yp,zp,t]
                            avgcount[x-p0[X],y-p0[Y],z-p0[Z],r] += 1
						
    # final step: average, discard data with unsufficient number of subjects
    for x in xrange(p0[X],pM[X]):
        for y in xrange(p0[Y],pM[Y]):
            for z in xrange(p0[Z],pM[Z]):
                for r in xrange(nruns) :
                    if (avgcount[x-p0[X],y-p0[Y],z-p0[Z],r]<=groupcutoff*nsubjects) :
                        avgmask[x,y,z] = 0
                
                if (avgmask[x,y,z]>0) :
                    for r in xrange(nruns) :
                        for t in xrange(runtimes[r], runtimes[r+1]) :
                            avgfmri[x-p0[X],y-p0[Y],z-p0[Z],t] /= avgcount[x-p0[X],y-p0[Y],z-p0[Z],r]
    if saveavgfmri :
        imgname = os.path.basename(fmrifiles[0])
        basename = imgname.split(".")
        outname = outdir+basename[0]+"_avgfmri."+'.'.join(basename[1:])
        print("saving to: "+outname)
        outimg = nb.Nifti1Image(avgfmri, affine=None)
        outimg.to_filename(outname)	
        
        outname = outdir+basename[0]+"_avgcount."+'.'.join(basename[1:])
        print("saving to: "+outname)
        outimg = nb.Nifti1Image(avgcount, affine=None)
        outimg.to_filename(outname)	
        
        outname = outdir+basename[0]+"_avgmask."+'.'.join(basename[1:])
        print("saving to: "+outname)
        outimg = nb.Nifti1Image(avgmask, affine=None)
        outimg.to_filename(outname)	
            
    return (avgfmri, avgcount, avgmask, p0, pM, runtimes)


def build_zscore(fmri, mask, transform, skip_first, skip_last, meancutoff) :
    nix = fmri.shape[X]
    niy = fmri.shape[Y]
    niz = fmri.shape[Z]
    # find mean, stdev,max within mask
    for xi in xrange(nix):
        for yi in xrange(niy):
            for zi in xrange(niz):
                fmean = 0
                fmax = 0
                fden = 0
                if (mask[xi,yi,zi]>0) :
                    for ti in xrange(skip_first, fmri.shape[T]-skip_last) :
                        fmean += fmri[xi,yi,zi,ti]
                        if (fmri[xi,yi,zi,ti]>fmax) : fmax = fmri[xi,yi,zi,ti]
                        fden += 1
                if (fden>0) : fmean /= fden
                
                # mask out data with low raw values
                if (fmean<=meancutoff*fmax) : mask[xi,yi,zi] = 0
                
                fstd = 0
                if (mask[xi,yi,zi]>0) :
                    for ti in xrange(skip_first, fmri.shape[T]-skip_last) :
                        fmri[xi,yi,zi,ti] -= fmean
                        fstd += fmri[xi,yi,zi,ti]*fmri[xi,yi,zi,ti]
                if (fden>1) : fstd = np.sqrt(fstd/(fden-1))
                
                # z-score
                if (mask[xi,yi,zi]>0) and (fstd>0) :
                    for ti in xrange(skip_first, fmri.shape[T]-skip_last) :
                        fmri[xi,yi,zi,ti] /= fstd
                else :
                    for ti in xrange(skip_first, fmri.shape[T]-skip_last) :
                        fmri[xi,yi,zi,ti] = 0
    return fmri

main()
