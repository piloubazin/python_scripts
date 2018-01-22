#!/usr/bin/python

# adapted from SBL bandpass filtering toolobox in Matlab

import nibabel as nb
import numpy as np
import sys
import os
import argparse

# global variables / flags

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

    parser = argparse.ArgumentParser()
    parser.add_argument('command')
    parser.add_argument('-o', '--outdir', nargs=1, help='output directory (default: input directory)')
    parser.add_argument('-f', '--fmri', nargs='*', metavar='sub#:run1,run2...', help='list of fMRI time series (4D) per run per subject')
    parser.add_argument('-m', '--mask', nargs='*', metavar='sub#:run1,run2...', help='list of masks (3D) per run per subject')
    parser.add_argument('-t', '--transform', nargs='*', metavar='sub#:run1,run2...', help='list of transformations (coordinate mappings) to group space per run per subject')
    parser.add_argument('-jf', '--jsonfile', nargs=1, help='json file listing fmri data, masks, and transforms')
    parser.add_argument('-mc', '--meancutoff', nargs='?', default=0.0001, help='mask out voxels with (normalized) mean over time below cutoff')
    parser.add_argument('-gc', '--groupcutoff', nargs='?', default=0.5, help='discard voxels with fewer subjects than this fraction of the group')
    parser.add_argument('-sf', '--skip_first', nargs='?', default=0, help='skip time frames at the beginning of each run')
    parser.add_argument('-sl', '--skip_last', nargs='?', default=0, help='skip time frames at the end of each run')
    parser.add_argument('-nr', '--nruns', nargs='?', default=1, help='number of runs per subject')
    print(parser.parse_args(sys.argv))
    
    # load the parameters
    subject_list = None
    json_file = parser.parse_args(sys.argv).jsonfile
    if json_file != None:
        inputs = input_json_file(json_file[0])
        fmri_files = inputs['fmri']
        mask_files = inputs['mask']
        transform_files = inputs['transform']
        subject_list = inputs['subject']
        run_list = inputs['run']
    else:
        fmri_files = parser.parse_args(sys.argv).fmri
        mask_files = parser.parse_args(sys.argv).mask
        transform_files = parser.parse_args(sys.argv).transform
        
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
    
    if subject_list != None :
        nsubjects = max(subject_list)
        nruns = max(run_list)
    else:    
        nruns = int(parser.parse_args(sys.argv).nruns)
        nsubjects = int(len(fmri_files)/nruns)
        run_list = []
        subject_list = []
        for sub in range(nsubjects):
            for run in range(nruns):
                run_list.append(int(run)+1)
                subject_list.append(int(sub)+1)
    
    print("subjects = "+str(subject_list))
    print("runs = "+str(run_list))
    
    inter_subject_correlation(fmri_files, mask_files, transform_files, outdir, nsubjects, subjects=subject_list,
                                nruns=nruns, runs=run_list, skip_first=skip_first, skip_last=skip_last, 
                                meancutoff=meancutoff, groupcutoff=groupcutoff)
    return


def inter_subject_correlation(fmri_files, mask_files, transform_files, outdir, nsubjects, subjects=None,
                                nruns=1, runs=None, skip_first=3, skip_last=0, meancutoff=0.0001, groupcutoff=0.5) :

    # 0. setup global variables and such
    basename = os.path.basename(fmri_files[0])
    basename = basename.split(".")

    # if subjects,runs not set build from numbers
    if (subjects==None) :         
        subjects = []
        for sub in range(nsubjects):
            for run in range(nruns):
                 subjects.append(int(sub)+1)

    if (runs==None) :         
        runs = []
        for sub in range(nsubjects):
            for run in range(nruns):
                runs.append(int(run)+1)

    print('ISC: '+str(nsubjects)+" subjects x "+str(nruns)+" runs, (mc,gc,sf,sl): ",meancutoff,groupcutoff,skip_first,skip_last)

    # 1. build a global average of all subjects (removing voxels with low intensity, missing values, and z-scoring the rest)
    (avgfmri, avgcount, avgmask, p0, pM, runtimes) = build_average(fmri_files, mask_files, transform_files, outdir,
                                                                    nsubjects, subjects, nruns, runs,
                                                                    skip_first, skip_last, meancutoff, groupcutoff)
    
    # 2. for each subject, remove data from the average and compute the correlation
    correlation = compute_subject_correlations(fmri_files, mask_files, transform_files, outdir,
                                               nsubjects, subjects, nruns, runs, 
                                               skip_first, skip_last, meancutoff,
                                               avgfmri, avgcount, avgmask, p0, pM, runtimes)
    
    # export results (resampled to avg space)
    transform = nb.load(transform_files[0])
    nx = transform.shape[X]
    ny = transform.shape[Y]
    nz = transform.shape[Z]
    resampled = np.zeros((nx,ny,nz,nsubjects))
    resampled[p0[X]:pM[X],p0[Y]:pM[Y],p0[Z]:pM[Z],:] = correlation
#    for x in range(p0[X],pM[X]):
#        for y in range(p0[Y],pM[Y]):
#            for z in range(p0[Z],pM[Z]):
#                if (avgmask[x,y,z]>0) :
#                    for s in range(nsubjects):
#                        resampled[x,y,z,s] = correlation[x-p0[X],y-p0[Y],z-p0[Z],s]
                        
    outname = outdir+basename[0]+"_isc"+str(nsubjects)+"x"+str(nruns)+"_corr."+".".join(basename[1:])
    print("saving to: "+outname)
    header = transform.get_header()
    affine = transform.get_affine()
    header['cal_min'] = np.min(resampled)
    header['cal_max'] = np.max(resampled)
    outimg = nb.Nifti1Image(resampled, affine, header)
    outimg.to_filename(outname)
    
    # 3. additional statistics?

    # outputs the filename rather than the entire file itself for easier handling
    return outname


def compute_subject_correlations(fmri_files, mask_files, transform_files, outdir,
                                nsubjects, subjects, nruns, runs,
                                skip_first, skip_last, meancutoff,
                                avgfmri, avgcount, avgmask, p0, pM, runtimes) :

    # test for consistency
    # here we assume the same numbers of subjects x runs for all data
    # (different versions could be made for simpler cases)
    if not ( len(fmri_files) == len(mask_files) 
            and len(fmri_files) == len(transform_files) 
            and len(fmri_files) == len(subjects) 
            and len(fmri_files) == len(runs) ):
        print("!different number of masks, transform and fmri data than specified subjects and runs!")
        return

    # for each subject, remove from average and correlate
    if debug : print("inter-subject correlations")    
    correlation = np.zeros((pM[X]-p0[X],pM[Y]-p0[Y],pM[Z]-p0[Z],nsubjects))
    for s in range(nsubjects) :
        if debug : print("subject: "+str(s+1))    

        # keep sums in case of multiple runs per subject
        
        covar = np.zeros((pM[X]-p0[X],pM[Y]-p0[Y],pM[Z]-p0[Z]), dtype=np.float64)
        var_sub = np.zeros((pM[X]-p0[X],pM[Y]-p0[Y],pM[Z]-p0[Z]), dtype=np.float64)
        var_avg = np.zeros((pM[X]-p0[X],pM[Y]-p0[Y],pM[Z]-p0[Z]), dtype=np.float64)
        samples = np.zeros((pM[X]-p0[X],pM[Y]-p0[Y],pM[Z]-p0[Z]), dtype=np.float64)

        for r in range(nruns) :
            if debug : print("run: "+str(r+1))    

            # search for (subject,run) combination
            idx = -1
            for n in range(len(fmri_files)) :
                if (subjects[n]==s+1 and runs[n]==r+1) : 
                    idx = n
                    break
            if (idx>-1) :
                if debug : print("data set: "+str(idx+1))    
                
                subjectrun = nb.load(fmri_files[idx])
                fmri = subjectrun.get_data()
                transform = nb.load(transform_files[idx]).get_data()
                mask = nb.load(mask_files[idx]).get_data()
            
                fmri = build_zscore(fmri, mask, transform, skip_first, skip_last, meancutoff)
                
                # remove from global average
#                for x in range(p0[X],pM[X]):
#                    for y in range(p0[Y],pM[Y]):
#                        for z in range(p0[Z],pM[Z]):
#                            if (avgmask[x,y,z]>0) :
#                                xp = int(np.rint(transform[x,y,z,X]))
#                                yp = int(np.rint(transform[x,y,z,Y]))
#                                zp = int(np.rint(transform[x,y,z,Z]))
#                            
#                                if (mask[xp,yp,zp]>0) : # only where you have valid data
#                                   for t in range(skip_first, fmri.shape[T]-skip_last) :
#                                       # remove normalization
#                                        avgdata = avgfmri[x-p0[X],y-p0[Y],z-p0[Z],runtimes[r]+t-skip_first] \
#                                                 *avgcount[x-p0[X],y-p0[Y],z-p0[Z],r]
#                                        # remove data value
#                                        avgdata -= fmri[xp,yp,zp,t]
#                                        # re-normalize
#                                        avgdata /= (avgcount[x-p0[X],y-p0[Y],z-p0[Z],r]-1)
#                                        # update variances and covariance
#                                        covar[x-p0[X],y-p0[Y],z-p0[Z]] += fmri[xp,yp,zp,t]*avgdata
#                                        var_sub[x-p0[X],y-p0[Y],z-p0[Z]] += fmri[xp,yp,zp,t]*fmri[xp,yp,zp,t]
#                                        var_avg[x-p0[X],y-p0[Y],z-p0[Z]] += avgdata*avgdata
#                                        samples[x-p0[X],y-p0[Y],z-p0[Z]] += 1

                # remove from global average
                transform = np.rint(transform).astype(int)
                cropped = transform[p0[X]:pM[X],p0[Y]:pM[Y],p0[Z]:pM[Z],:]
                # remove current subject from average
                avgdata = avgfmri[:,:,:,runtimes[r]:runtimes[r+1]]*avgcount[:,:,:,r,np.newaxis]
                avgdata -= fmri[cropped[:,:,:,X],cropped[:,:,:,Y],cropped[:,:,:,Z],skip_first:fmri.shape[T]-skip_last]
                avgdata /= np.fmax(avgcount[:,:,:,r,np.newaxis]-1,np.ones(avgdata.shape))
                # update variance and covariance
                for t in range(skip_first, fmri.shape[T]-skip_last) :
                    covar += np.multiply(fmri[cropped[:,:,:,X],cropped[:,:,:,Y],cropped[:,:,:,Z],t],avgdata[:,:,:,t-skip_first])
                    var_sub += np.multiply(fmri[cropped[:,:,:,X],cropped[:,:,:,Y],cropped[:,:,:,Z],t],fmri[cropped[:,:,:,X],cropped[:,:,:,Y],cropped[:,:,:,Z],t])
                    var_avg += np.multiply(avgdata[:,:,:,t-skip_first],avgdata[:,:,:,t-skip_first])    
                    samples += 1

        # compute the final correlations
#        for x in range(p0[X],pM[X]):
#            for y in range(p0[Y],pM[Y]):
#                for z in range(p0[Z],pM[Z]):
#                    if (avgmask[x,y,z]>0) :
#                        correlation[x-p0[X],y-p0[Y],z-p0[Z],s] = covar[x-p0[X],y-p0[Y],z-p0[Z]] \
#                                                                /np.sqrt(var_sub[x-p0[X],y-p0[Y],z-p0[Z]]) \
#                                                                /np.sqrt(var_avg[x-p0[X],y-p0[Y],z-p0[Z]])
 
        # compute the final correlations
        correlation[:,:,:,s] = covar/(np.sqrt(np.multiply(var_sub,var_avg)))
 
    return correlation                       


def build_average( fmri_files, mask_files, transform_files, outdir,
                    nsubjects, subjects, nruns, runs, skip_first, skip_last, 
                    meancutoff, groupcutoff) :
    
    # test for consistency
    # here we assume the same numbers of subjects x runs for all data
    # (different versions could be made for simpler cases)
    if not ( len(fmri_files) == len(mask_files) 
            and len(fmri_files) == len(transform_files) 
            and len(fmri_files) == len(subjects) 
            and len(fmri_files) == len(runs) ):
        print("!different number of masks, transform and fmri data than specified subjects and runs!")
        return
        
    # init from first subject
    if debug : print("opening file: "+mask_files[0])
    mask = nb.load(mask_files[0]).get_data()
    
    if debug : print("opening file: "+transform_files[0])
    transform = nb.load(transform_files[0]).get_data()
    
    # define common space: every non-zero element of the mask is used
    nx = transform.shape[X]
    ny = transform.shape[Y]
    nz = transform.shape[Z]
    
    avgmask = np.zeros(transform.shape[X:T])
    for n in range(len(mask_files)) :
        if debug : print("subject: "+str(subjects[n])+", run:"+str(runs[n]))
        transform = nb.load(transform_files[n]).get_data()
        mask = nb.load(mask_files[n]).get_data()

        transform = np.rint(transform).astype(int)
        avgmask += mask[transform[:,:,:,X],transform[:,:,:,Y],transform[:,:,:,Z]]
             
    # find mask boundaries
    indices = np.nonzero(avgmask)
    p0 = np.array([int(np.min(indices[X])),int(np.min(indices[Y])),int(np.min(indices[Z]))])
    pM = np.array([int(np.max(indices[X])+1),int(np.max(indices[Y])+1),int(np.max(indices[Z])+1)])
    print("analysis bounding box (global space): "+str(p0)+", "+str(pM))
    
    # create average fmri only within the mask, concatenate the runs
    runtimes = list(range(nruns+1))
    # first run time is zero, simpler for later computations
    runtimes[0] = 0
    for n in range(len(fmri_files)) :
        f = runs[n]
        fmri_time = nb.load(fmri_files[n]).shape[T]
        runtimes[f] = int(runtimes[f-1] + fmri_time-skip_first-skip_last)
        
    if debug : print("runtimes: "+str(runtimes))    
            
    avgfmri = np.zeros((pM[X]-p0[X],pM[Y]-p0[Y],pM[Z]-p0[Z],runtimes[nruns]))
    # to count the number of subjects per run that contribute to each voxel
    avgcount = np.zeros((pM[X]-p0[X],pM[Y]-p0[Y],pM[Z]-p0[Z],nruns))
        
    # pull data from each time series
    # check for % missing data (below threshold) and mask out the bad ones
    # use the rest to build Z-scores
    # note: for simplicity, the number of missing data points is the maximum over runs per voxel
    if debug : print("global averaging")    
    for n in range(len(fmri_files)) :
        s = subjects[n]
        r = runs[n]

        if debug : print("subject: "+str(subjects[n])+", run: "+str(runs[n]))    

        subjectrun = nb.load(fmri_files[n])
        fmri = subjectrun.get_data()
        transform = nb.load(transform_files[n]).get_data()
        mask = nb.load(mask_files[n]).get_data()
    
        fmri = build_zscore(fmri, mask, transform, skip_first, skip_last, meancutoff)
        # option: save the masked, z-scored time series?
        if savezscores :
            imgname = os.path.basename(fmri_files[n])
            basename = imgname.split(".")
            outname = outdir+basename[0]+"_zscored."+'.'.join(basename[1:])
            print("saving to: "+outname)
            header = subjectrun.get_header()
            affine = subjectrun.get_affine()
            header['cal_min'] = np.min(fmri)
            header['cal_max'] = np.max(fmri)
            outimg = nb.Nifti1Image(fmri, affine, header)
            outimg.to_filename(outname)	

        # build the global average
#        for x in range(p0[X],pM[X]):
#            for y in range(p0[Y],pM[Y]):
#                for z in range(p0[Z],pM[Z]):
#                    xp = int(np.rint(transform[x,y,z,X]))
#                    yp = int(np.rint(transform[x,y,z,Y]))
#                    zp = int(np.rint(transform[x,y,z,Z]))
#                    
#                    if (mask[xp,yp,zp]>0) :
#                        for t in range(skip_first, fmri.shape[T]-skip_last) :
#                            avgfmri[x-p0[X],y-p0[Y],z-p0[Z],runtimes[r-1]+t-skip_first] += fmri[xp,yp,zp,t]
#                        avgcount[x-p0[X],y-p0[Y],z-p0[Z],r-1] += 1
						
        # build the global average: array version?
        transform = np.rint(transform).astype(int)
        cropped = transform[p0[X]:pM[X],p0[Y]:pM[Y],p0[Z]:pM[Z],:]
        #np.clip(cropped[:,:,:,X], 0, pM[X]-p0[X]-1, out=cropped[:,:,:,X])
        #np.clip(cropped[:,:,:,Y], 0, pM[Y]-p0[Y]-1, out=cropped[:,:,:,Y])
        #np.clip(cropped[:,:,:,Z], 0, pM[Z]-p0[Z]-1, out=cropped[:,:,:,Z])
        
        avgfmri[:,:,:,runtimes[r-1]:runtimes[r]] \
                += (mask[cropped[:,:,:,X],cropped[:,:,:,Y],cropped[:,:,:,Z],np.newaxis]>0) \
                *fmri[cropped[:,:,:,X],cropped[:,:,:,Y],cropped[:,:,:,Z],skip_first:fmri.shape[T]-skip_last]
                
        avgcount[:,:,:,r-1] \
                += (mask[cropped[:,:,:,X],cropped[:,:,:,Y],cropped[:,:,:,Z]]>0)
						
    # final step: average, discard data with unsufficient number of subjects
#    for x in range(p0[X],pM[X]):
#        for y in range(p0[Y],pM[Y]):
#            for z in range(p0[Z],pM[Z]):
#                for r in range(nruns) :
#                    if (avgcount[x-p0[X],y-p0[Y],z-p0[Z],r]<=groupcutoff*nsubjects) :
#                        avgmask[x,y,z] = 0
#                
#                if (avgmask[x,y,z]>0) :
#                    for r in range(nruns) :
#                        for t in range(runtimes[r], runtimes[r+1]) :
#                            avgfmri[x-p0[X],y-p0[Y],z-p0[Z],t] /= avgcount[x-p0[X],y-p0[Y],z-p0[Z],r]

    # final step: average, discard data with unsufficient number of subjects
    for r in range(nruns) :
        avgmask[p0[X]:pM[X],p0[Y]:pM[Y],p0[Z]:pM[Z]] *= (avgcount[:,:,:,r]>groupcutoff*nsubjects)

    np.fmax(avgcount,np.ones(avgcount.shape), out=avgcount)
    for r in range(nruns) :
        avgfmri[:,:,:,runtimes[r]:runtimes[r+1]] /= avgcount[:,:,:,r,np.newaxis]

    if saveavgfmri :
        imgname = os.path.basename(fmri_files[0])
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
#    nix = fmri.shape[X]
#    niy = fmri.shape[Y]
#    niz = fmri.shape[Z]
#    # find mean, stdev,max within mask
#    for xi in range(nix):
#        for yi in range(niy):
#            for zi in range(niz):
#                fmean = 0
#                fmax = 0
#                fden = 0
#                if (mask[xi,yi,zi]>0) :
#                    for ti in range(skip_first, fmri.shape[T]-skip_last) :
#                        fmean += fmri[xi,yi,zi,ti]
#                        if (fmri[xi,yi,zi,ti]>fmax) : fmax = fmri[xi,yi,zi,ti]
#                        fden += 1
#                if (fden>0) : fmean /= fden
#                
#                # mask out data with low raw values
#                if (fmean<=meancutoff*fmax) : mask[xi,yi,zi] = 0
#                
#                fstd = 0
#                if (mask[xi,yi,zi]>0) :
#                    for ti in range(skip_first, fmri.shape[T]-skip_last) :
#                       fmri[xi,yi,zi,ti] -= fmean
#                       fstd += fmri[xi,yi,zi,ti]*fmri[xi,yi,zi,ti]
#                if (fden>1) : fstd = np.sqrt(fstd/(fden-1))
#                
#                # z-score
#                if (mask[xi,yi,zi]>0) and (fstd>0) :
#                    for ti in range(skip_first, fmri.shape[T]-skip_last) :
#                        fmri[xi,yi,zi,ti] /= fstd
#                else :
#                    for ti in range(skip_first, fmri.shape[T]-skip_last) :
#                        fmri[xi,yi,zi,ti] = 0

    # find mean, stdev,max within mask
    fmean = np.mean(fmri[:,:,:,skip_first:fmri.shape[T]-skip_last], axis=T, dtype=np.float64)
    fstd = np.std(fmri[:,:,:,skip_first:fmri.shape[T]-skip_last], axis=T, dtype=np.float64)
    fmax = np.max(fmean)
    # mask out data with low raw values
    mask *= (fmean>meancutoff*fmax)      
    # z-score
    fmri -= fmean[:,:,:,np.newaxis]
    fmri = np.divide(fmri, fstd[:,:,:,np.newaxis], where=(mask>0)[:,:,:,np.newaxis])
    fmri *= (mask>0)[:,:,:,np.newaxis]

    return fmri


def input_json_file(json_file) :
    import json
    
    print("Loading :"+json_file)
    
    with open(json_file, 'r') as f:
        data = json.load(f)

    basedir = os.path.dirname(json_file)

    subject = []
    run = []
    fmri = []
    mask = []
    transform = []
    for record in data :
        subject.append(int(record['subject']))
        run.append(int(record['run']))
        fmri.append(os.path.join(basedir, str(record['fmri'])))
        mask.append(os.path.join(basedir, str(record['mask'])))
        transform.append(os.path.join(basedir, str(record['transform'])))
 
    #print("subjects = "+str(subject))
    #print("runs = "+str(run))
    #print("fmri = "+str(fmri))
    #print("mask = "+str(mask))
    #print("transform = "+str(transform))

    # test for existence of the files
    for path in fmri :
        if not os.path.isfile(path) :
            print("missing file: "+str(path))
            return
            
    for path in mask :
        if not os.path.isfile(path) :
            print("missing file: "+str(path))
            return
            
    for path in transform :
        if not os.path.isfile(path) :
            print("missing file: "+str(path))
            return
    
    print("all files accounted for")
    
    return {'fmri':fmri, 'mask':mask, 'transform':transform, 'subject':subject, 'run':run}


if __name__ == "__main__":
    main()