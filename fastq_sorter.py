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
        
def run_cmd(lis):
    print(" ".join(lis))
    os.system(" ".join(lis))

outfile='/home/javiernunezgarcia/files_to_delete.csv'
outfile_summary='/home/javiernunezgarcia/files_to_delete_summary.csv'
mother_paths=['/home/javiernunezgarcia/mnt/fsx-016',\
              '/home/javiernunezgarcia/mnt/fsx-024',\
              '/home/javiernunezgarcia/mnt/fsx-ranch-023']

fils_ends=['_alignment_stats.csv.gz',\
           '.pileup.vcf.gz',\
           '.pileup.vcf',\
           '.flt.snp', \
           '.snp',\
           '.sorted.bam',\
           '.sorted.bam.bai',\
           '_DeletedSequences.csv',\
           '_R1_fastqc.zip',\
           '_R2_fastqc.zip',\
           'Mapped_fastqs_Info.txt']  

#'_alignment_stats.csv'

tab=[["type","size in MB","file"]]
s_tab=[["type","path","size in MB","#files"]]

for fil_end in fils_ends:
    for mother_path in mother_paths:
        fils=glob.glob(mother_path+'/**/*'+fil_end, recursive=True)
        sub_tab=[]
        for fil in fils:
            sub_tab.append([fil_end,round(os.path.getsize(fil)/float((1024*1024)),2),fil])
        print(fil_end+"  "+mother_path+"   "+str(sum([x[1] for x in sub_tab]))+"   "+str(len(fils)))
        s_tab.append([fil_end,mother_path,sum([x[1] for x in sub_tab]),str(len(fils))])    
        tab=tab+sub_tab

writecsv(outfile,tab)
writecsv(outfile_summary,s_tab)
    


