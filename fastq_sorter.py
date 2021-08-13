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
    mother_path='/home/javiernunezgarcia/mnt/fsx-016'

#_alignment_stats.csv.gz
#.pileup.vcf.gz
#.flt.snp
#.sorted.bam
#.sorted.bam.bai    
    
#### making a list of all the fastq files (faster if we look each sample one by one)
pileup_fils=glob.glob(mother_path+'/**/*.pileup.vcf.gz', recursive=True)
alignment_fils=glob.glob(mother_path+'/**/*_alignment_stats.csv.gz', recursive=True)
fltsnp_fils=glob.glob(mother_path+'/**/*.flt.snp', recursive=True)
sortedbam_fils=glob.glob(mother_path+'/**/*.sorted.bam', recursive=True)
sortedbambai_fils=glob.glob(mother_path+'/**/*.sorted.bam.bai', recursive=True)



