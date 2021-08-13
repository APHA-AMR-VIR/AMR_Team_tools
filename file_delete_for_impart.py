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
samples_tab=[x[3] for x in samples_tab[1:] if x[3]!=""]
print(len(samples_tab))
results=[]
cont=0
for fil in samples_tab:
    print(cont)
    cont=cont+1
    r1=fil
    if "_R1_" in r1:
        r2=fil.replace("_R1_","_R2_")
        run_cmd(["rm","\""+r1+"\""])
        run_cmd(["rm","\""+r2+"\""])
        if not os.path.isfile(r1) and not os.path.isfile(r2):
            results.append([r1,r2,"deleted"])
        else:
            results.append([r1,r2,"still_alive"])
    else:
        print("***************************************************")
        print(r1)
        print("****************************************************")
        results.append([r1,r2,"no_pattern"])

results=[["r1","r2","status"]]+results       
writecsv(os.path.join('/home/javiernunezgarcia','summary_deletion.csv'),results)