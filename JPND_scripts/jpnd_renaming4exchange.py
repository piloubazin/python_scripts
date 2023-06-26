# simple script to select specific data form a BIDS-like database, 
# take a random subset, and rename for transfer to partner institution.
#
# Important:
# 1. All the data should be anonymized first (removal of identifying subject information, 
#    skull strip for MRI, etc). This should be done anyway when building the database.
#
# 2. This script expects a BIDS-like structure, which means every subject is in its own directory
#    and all corresponding data to search is within that directory. The way things are organized
#    within each subject folder is not important.
#
# The script will identify all subjects with the desired modalities, randomize subject order,
# build a new folder with new subject naming but identical BIDS-like structure, 
# copy the data over. 
# The randomization key is *not* saved anywhere, so the data is fully GDPR compliant.

import os
import glob
import sys
import random
import datetime
import shutil

# short names for each site (do not change!)
sites = {"leipzig": "lz", "amsterdam": "am", "liege": "li", "pecs": "pe"}

# site-specific parameters:

# source center (fixed)
source = sites["liege"]

# location of the database on disk (fixed)
data_dir = "/home/Public/jpnd/data/Liege@UVA/scratch/ulg/crc/ghammad/GWA/EPI/"


# target center (changing)
target = sites["amsterdam"]

# number of subjects to transfer (changing)
ndata = 6

# type of image modalities to transfer (changing)
# (each name has to be uniquely descriptive of a single file per subject)
modalities = ["R1",\
              "R2s_OLS",\
              "PD",\
			  "MTsat"]
              
# data type (changing)
data_types = [".nii", \
              ".nii", \
              ".nii", \
              ".nii"]

# corresponding sub-folder structure (changing)
# (under the main 
subfolders = ["AutoReorient/Results/",\
              "AutoReorient/Results/",\
              "AutoReorient/Results/",\
              "AutoReorient/Results/"]
           
# where to write the resulting data set
output_dir = "/home/Public/jpnd/data/liege/"


# beyond this everything should remain the same (do not change!)

# step 1: look for subject subfolders
subjects = sorted(glob.glob(data_dir+"*/"))
if len(subjects)==0:
    print("no data found, check path: "+data_dir)
    sys.exit()
    
print("found "+str(len(subjects))+" subject subfolders")


# step 2: build a list of all subjects
for idx,subject in enumerate(subjects):
    subject = subject.replace(data_dir,"")
    subjects[idx] = subject.replace("/","")

# step 3: look for desired data, keep only subjects with all modalities
for subject in subjects:
    for idx,modality in enumerate(modalities):
        search = sorted(glob.glob(os.path.join(data_dir,subject,subfolders[idx],"*"+modality+"*"+data_types[idx])))

        if len(search)==0:
            print("no data found for subject ("+os.path.join(data_dir,subject,subfolders[idx],"*"+modality+"*"+data_types[idx])+"), remove from candidates")
            if subjects.count(subject)>0: subjects.remove(subject)
        if len(search)>1:
            print("multiple data found for subject ("+os.path.join(data_dir,subject,subfolders[idx],"*"+modality+"*"+data_types[idx])+"), remove from candidates")
            if subjects.count(subject)>0: subjects.remove(subject)
            
print("subjects with all required data: "+str(len(subjects)))

if (len(subjects)<=ndata):
    print("not enough subjects for randomization!")
    sys.exit()


# step 4: shuffle remaining data, pick desired number, build new file subfolders
indices = list(range(len(subjects)))
random.shuffle(indices)

today = datetime.datetime.now().strftime("%y%m%d%H%M%S")

for num,sub in enumerate(indices[0:ndata]):
    # create new subject name with sub-source2targetdatennum as subject id
    newsub = 'sub-'+source+'2'+target+today+"n"+str(num+1).zfill(4)
    
    # build main directory
    os.makedirs(os.path.join(output_dir,newsub))
    
    for idx,modality in enumerate(modalities):
        # copy full source directory structure
        os.makedirs(os.path.join(output_dir,newsub,subfolders[idx]), exist_ok=True)
        # copy data to new name
        data = glob.glob(os.path.join(data_dir,subjects[sub],subfolders[idx],"*"+modality+"*"+data_types[idx]))[0]
        copied = os.path.join(output_dir,newsub,subfolders[idx],newsub+"_"+modality+data_types[idx])
        copied = copied.replace(subjects[sub], newsub)
        print("write "+os.path.join(newsub,subfolders[idx],newsub+"_"+modality+data_types[idx]))
        shutil.copyfile(data, copied)
    
print("done")    

