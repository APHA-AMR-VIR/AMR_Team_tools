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

args=sys.argv
if len(args)>1:
    mother_path=args[1]
    outfile=args[2]
else:
    outfile='/home/javiernunezgarcia/fastqs_at_IMPART_2018.csv'
    mother_path='/home/javiernunezgarcia/mnt/fsx-016/Non-Defra_Projects/EU_Projects/IMPART_2018'

#### making a list of all the fastq files (faster if we look each sample one by one)
fils=glob.glob(mother_path+'/**/*R1*.fastq.gz', recursive=True)
samples=set([fil.split(os.sep)[-1].split("_")[0] for fil in fils])

tab=[["sample","R1 fastqs"]]
for sam in samples:
    tab.append([sam,""])
    sub_fils=[fil for fil in fils if fil.split(os.sep)[-1].split("_")[0]==sam]
    for s in sub_fils:
        tab.append(["",s,int(os.path.getsize(s))])
    
writecsv(outfile,tab)             


