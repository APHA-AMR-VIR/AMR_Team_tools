# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 10:53:42 2021

@author: Javier Nunez


"""

import os,os.path,sys
import glob,csv


####### functions definitions
def readTable(fname,ch=','):
    infile=open(fname,"r")
    data = csv.reader(infile, delimiter=ch)
    dataOut = [row for row in data]
    infile.close()
    return(dataOut)

def writecsv(fname,matrix):
    with open(fname, "w") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)
        print("file "+fname+" saved.")
        
def run_cmd(lis):   ### this funtion takes a list of strings, concatenate them and run a command line
    print(" ".join(lis))
    os.system(" ".join(lis))


#### reading the arguments
   
#samples_tab_file=sys.argv[1]  ## table with the samples to be found
#mother_path=sys.argv[2]  ## path where to look for the fastq files
#folder_out=sys.argv[3]   ## folder where to copy the files


samples_tab_file='/home/javiernunezgarcia/samples_to_upload.csv'
mother_path='/home/javiernunezgarcia/mnt/fsx-016' #/RDVM0529/VM0529'
folder_out='~/upload_try'

#### making a list of all the fastq files (faster if we look each sample one by one)
fastqs=glob.glob(mother_path+'/**/*.fastq.gz', recursive=True)
print(str(len(fastqs))+" fastq files in path "+mother_path)

#### reading the table intoa list of lists (every line in the table is a list of items)
#### better option is using pandas!
samples_tab=readTable(samples_tab_file)


### creating the output folder
if not os.path.exists(folder_out):
    run_cmd(['mkdir','-p',folder_out])


### initialise a new table (with the columns names) that will contain the original table plus the fastq files location
### again pandas is a better option to manage tables
new_samples_tab=[samples_tab[0]+["number_R1","R1 sizes","number_R2","R2_sizes","R1s","R2s"]]

#### loop every line in the table to find the corresponding R1 and R2 files
#### the locations ar recorded in the new table and the files copied in the new destination folder

for row in samples_tab[1:]:
    R1=[x for x in fastqs if "_R1" in x and row[0]==x.split(os.sep)[-1][:len(row[0])]]
    R1_file_size=[round(os.path.getsize(x)/(1024*1024),4) for x in R1]
    R2=[x for x in fastqs if "_R2" in x and row[0]==x.split(os.sep)[-1][:len(row[0])]]
    R2_file_size=[round(os.path.getsize(x)/(1024*1024),4) for x in R2]
    new_samples_tab.append(row+[len(R1),R1_file_size,len(R2),R2_file_size,R1,R2])
    #if len(R1)==1 and len(R2)==1:
    #    run_cmd(["cp",R1[0],folder_out])
    #    run_cmd(["cp",R2[0],folder_out])
        
    
writecsv(samples_tab_file[:-4]+'_completed.csv',new_samples_tab)    



