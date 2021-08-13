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

samples_tab_file='/home/javiernunezgarcia/samples in SCE_v3.csv'
folder_out='/home/javiernunezgarcia/mnt/s3-ranch-022/WGS_data/Non-Defra_projects/EU_projects/IMPART'

samples_tab=readTable(samples_tab_file)
samples_tab=[x for x in samples_tab if x[-3]=="keep"]
print(len(samples_tab))
results=[]
cont=0
for row in samples_tab:
    print(cont)
    cont=cont+1
    r1=row[3]
    if "_R1_" in r1:
        r2=row[3].replace("_R1_","_R2_")
        r1_size=int(os.path.getsize(r1))
        r2_size=int(os.path.getsize(r2))
        spe=row[-1]
        dest_folder=os.path.join(folder_out,spe)
        if not os.path.exists(dest_folder):
            run_cmd(['mkdir','-p',dest_folder])
        r1_c=os.path.join(dest_folder,r1.split(os.sep)[-1])
        r2_c=os.path.join(dest_folder,r2.split(os.sep)[-1])
        if not os.path.isfile(r1_c) or not int(os.path.getsize(r1_c))==r1_size:
            run_cmd(["cp","\""+r1+"\"",r1_c])
        if not os.path.isfile(r2_c) or not int(os.path.getsize(r2_c))==r2_size:
            run_cmd(["cp","\""+r2+"\"",r2_c])
        r1_c_size=int(os.path.getsize(r1_c))
        r2_c_size=int(os.path.getsize(r2_c))
        if r1_size==r1_c_size:
            r1_copy="ok"
        else:
            r1_copy="fail"
        if r2_size==r2_c_size:
            r2_copy="ok"
        else:
            r2_copy="fail"
        results.append([r1,r1_c,r1_size,r1_c_size,r1_copy,r2,r2_c,r2_size,r2_c_size,r2_copy])
    else:
        print("***************************************************")
        print(r1)
        print("****************************************************")
        results.append([r1,"no_pattern","","","",r2,"no_pattern","","",""])

results=[["r1","r1_c","r1_size","r1_c_size","r1_copy","r2","r2_c","r2_size","r2_c_size","r2_copy"]]+results       
writecsv(os.path.join('/home/javiernunezgarcia','summary.csv'),results)