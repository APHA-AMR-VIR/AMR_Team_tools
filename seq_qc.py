"""
Created on Tue Mar 23 10:53:42 2021

@author: Javier Nunez

#Fast QC basic metrics for R1, R2 and contigs files
#Run:   python seq_qc.py /path/to/arguments/file/arguments_file.args

#The number of reads and depth of coverage should be as high as possible. 
#   There is no assessed cut-off numbers for these. The higher number of reads and depth of coverage 
#   indicate high amount of raw reads data to start with. The depth of coverage (C) can be calculated as the 
#   length of the reads (L) divided by the genome size (G) multiplied by the number of reads (N); C = N * (L/G). 
#The average read length should be similar to the expected read length from the selected sequencing platform. 
#The size of assembled genome should not deviate more than 0.5 million base-pairs from the expected genome size. 
#   If the deviation in the size of the assembled genome is greater than 0.5 million based-pair of expected genome 
#   size, this is an indication that the genome sequences are either contaminated, not the expected species or poor 
#   sequencing quality. 
#The total number of contigs (after assembly) should be less than 500 contigs. A higher number of contigs 
#   indicates poor sequencing quality. 
#N50 indicates size of contigs in general. The higher N50 indicates the longer contigs in a genome. 
#   There is no general cut-off for N50, but some suggest using a N50 of >30 000 bp (Bortolaia et al. 2020). 


"""
## install quast (conda or pip)
import gzip
#conda install -c conda-forge biopython
import os,os.path,fnmatch,csv,sys
from os import listdir
import numpy as np
from multiprocessing import Pool

def writeCSV(fname,matrix):
    with open(fname, "w") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)
        print("file "+fname+" saved.")

def getCMD(lis):
    print(" ".join(lis))
    os.system(" ".join(lis))
    
    
def find_file(pattern, path):
    result = []
    for f in os.listdir(path):
        if fnmatch.fnmatch(f, pattern):
            result.append(f) #os.path.join(root, name))
    return result

def qc_fastq(fname,len_ref):
    with gzip.open(fname, 'r') as fil:
        cont=1
        cont_seqs=0
        lengths=[]
        for line in fil:
            if cont%2==0 and cont%4!=0:
                lengths.append(len(line))
                cont_seqs=cont_seqs+1
            cont=cont+1
    return([cont_seqs,round(np.average(lengths),2)])

def fasta_len(fil):
    with open(fil, 'r') as fil:
        cont=1
        len_ref=0
        for line in fil:
            if cont!=1:
                len_ref=len_ref+len(line.strip())
            cont=cont+1
    return(len_ref)

def read_contigs(fname):
    fileIn = open(fname, 'r')
    lines= fileIn.readlines()
    lines=[x.strip() for x in lines]
    ids=[]
    seqs=[]
    seq=""
    for line in lines:
        if line[0]==">":
            if len(ids)>0:
                seqs.append(seq)
                seq=""
            ids.append(line[1:])
        else:
            seq=seq+line
    seqs.append(seq)
    return(ids,seqs)

def contigs_checker(fname,len_ref):
    ids,seqs=read_contigs(fname)
    lens=[len(seqs[i]) for i in range(len(seqs))]
    smaller_500=len([x for x in lens if x<500])
    assembly_len=sum(lens)
    ### N50
    cont=0
    sum_n50=0
    while sum_n50<=assembly_len/2:
        sum_n50=sum_n50+lens[cont]
        cont=cont+1
    n50=lens[cont-1]
    ### NG50
    cont=0
    sum_ng50=0
    while sum_ng50<=len_ref/2  and cont<len(lens):
        sum_ng50=sum_ng50+lens[cont]
        cont=cont+1
    ng50=lens[cont-1]
    return([len(seqs),smaller_500,assembly_len,min(lens),max(lens),round(np.mean(lens)),n50,ng50])

def filter_contigs(fin,fout,min_size=300):
    ids,seqs=read_contigs(fin)
    lines=[]
    for i in range(len(ids)):
        if len(seqs[i])>=300:
            if i==0:
                lines.append(">"+ids[i])
            else:
                lines.append("\n"+">"+ids[i])
            for j in range(0,len(seqs[i]),60):
                lines.append("\n"+seqs[i][j:j+60])

    fileOut = open(fout, 'w')
    fileOut.writelines(lines)
    fileOut.close()        


def one_sample(row):
    try:
        sample_name=row[0].split("_")[0]
        r1=os.path.join(fastqs_path,row[0])
        r2=os.path.join(fastqs_path,row[1])
        fasta=os.path.join(fastas_folder,row[2])
        contigs_res=contigs_checker(fasta,len_ref)
        if int(min_len)>0:
            filter_contigs(fasta,os.path.join(fastas_folder,row[2].split(".")[0]+"_filtered.fasta"),min_len)
        fastq_r1_res=qc_fastq(r1,len_ref)
        fastq_r2_res=qc_fastq(r2,len_ref)
        st_cov=round((fastq_r1_res[0]+fastq_r1_res[0])*np.mean([fastq_r1_res[1],fastq_r1_res[1]])/len_ref,2)
        return([sample_name,ref_name,len_ref,fasta.split(os.sep)[-1]]+contigs_res+[r1.split(os.sep)[-1]]+fastq_r1_res+[r2.split(os.sep)[-1]]+fastq_r2_res+[st_cov])

    except:
        print("Something went wrong.")
        print("Therefore, "+sample_name+" no processed at all.")
        return(19*["NA"])               


###########Command line arguments parshing
fastqs_path=""
fastas_folder=""
R1_pattern=""
reference_genome=""
out_file=""
ncores=1

args=sys.argv
if len(args)>1:
    arguments_file=args[1]

else:
    arguments_file="/home/javiernunezgarcia/AMR_Team_tools/seq_qc_arguments_template.args"
    print("Argument file not given or doesn't exist. Please re run with /full/path/to/arguments/file/arguments_file.args")
    #sys.exit()

########Loading other arguments from the arguments file
print('Reading arguments from file: '+arguments_file)
with open(arguments_file,'r') as f:
    for line in f:
        if line[0] not in ['\n',' ','#']:
            print('Loading argument: '+ line.strip())
            exec(line.strip())
  

print("***** Checking samples to be run")
fastq_R1s=find_file("*"+R1_pattern+"*.fastq.gz", fastqs_path)
R2_pattern=R1_pattern.replace("R1","R2")
fastq_to_process=[]

for fastq_R1 in fastq_R1s:
    fastq_R2=fastq_R1.replace(R1_pattern,R2_pattern)
    sample_name=fastq_R1.split(os.sep)[-1].split("_")[0]
    if os.path.isfile(os.path.join(fastqs_path,fastq_R2)):
        R2_ok="yes"
    else:
        R2_ok="no"
        print("Missing R2 fastq for file: "+fastq_R1)
    fasta_file=[f for f in listdir(fastas_folder) if f[:len(sample_name)]==sample_name and f.split(".")[-1] in ["fasta","fa"]]    
    if len(fasta_file)==1:
        fasta_ok="yes"
        fasta_file=fasta_file[0]
    elif len(fasta_file)==0:
        fasta_ok="no"
        print("Missing assembly file for file: "+fastq_R1)
    else:
        fasta_ok="no"
        print("More than one assembly for: "+fastq_R1)
       
    if R2_ok=="yes" and fasta_ok=="yes":
        fastq_to_process.append([fastq_R1,fastq_R2,fasta_file])
        
print("***** Processing "+str(len(fastq_to_process))+" samples.")


### initiating the output table 
tab=[["sample_name","reference","reference_length","Contigs_file","#contigs","#contigs<500","assembly_length","minimum_length","max_length","average_length","N50","NG50","R1_fastq","R1_#reads","R1_Average_read_length","R2_fastq","R2_#reads","R2_Average_read_length","Estimated_depth"]]

ref_name=reference_genome.split(os.sep)[-1] 
len_ref=fasta_len(reference_genome)
  
#for row in fastq_to_process:
#    one_sample(row)
    
pool=Pool(ncores)
result=pool.map(one_sample,fastq_to_process)

writeCSV(out_file,tab+result)













