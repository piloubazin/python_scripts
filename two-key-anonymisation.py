from pandas import DataFrame, read_csv
import os
import shutil
import glob

# load the key files
file1 = '/home/jalkema1/Ahead/Recode1.csv'
file2 = '/home/jalkema1/Ahead/Recode2.csv'
#file1 = '/home/pilou/Projects/tmp/test_keys1.csv'
#file2 = '/home/pilou/Projects/tmp/test_keys2.csv'

in_folder = '/home/public/ahead_1.0/'
out_folder = '/home/public/ahead_1.0_final/'
#in_folder = '/home/pilou/Projects/tmp/db'
#out_folder = '/home/pilou/Projects/tmp/newdb'

#key1 = read_csv(file1, header=None)
#key2 = read_csv(file2, header=None)
key1 = read_csv(file1)
key2 = read_csv(file2)

os.makedirs(out_folder,exist_ok=True)

for subject in range(key1.shape[0]):
    # find the match
    for match in range(key2.shape[0]):
        if key1.values[subject][1]==key2.values[match][0]:
            old_name = key1.values[subject][0]
            new_name = key2.values[match][1]
        #if key1.values[subject]==key2.index[match]:
            #old_name = key1.index[subject]
            #new_name = key2.values[match]
            
    # copy directory and content
    print(old_name+" -> "+new_name)
    
    data_files = glob.glob(in_folder+old_name+'/**/'+old_name+'*',recursive=True)
    for data in data_files:
        print(data+" -> "+data.replace(in_folder,out_folder).replace(old_name,new_name))
        #print("create dir: "+os.path.dirname(data.replace(in_folder,out_folder).replace(old_name,new_name)))
        os.makedirs(os.path.dirname(data.replace(in_folder,out_folder).replace(old_name,new_name)),exist_ok=True)
        shutil.copyfile(data,data.replace(in_folder,out_folder).replace(old_name,new_name))
    
