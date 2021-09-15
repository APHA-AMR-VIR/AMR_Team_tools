### ASk Javier Nunez for any question
## This script runs unicyler in parallel. It uses the instance hard drive as temporary/auxiliar storage and 
## then it copy the assemblies where they come from in a folder caller "uc_fastas".
## Check the arguments file associted (same name) for info about the arguments you need to pass
 

import sys,os,os.path,fnmatch,csv
import multiprocessing as mp
from os import listdir
import random
from datetime import datetime
from pathlib import Path



def writeCSV(fname,matrix):
    with open(fname, "w") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)
        print("file "+fname+" saved.")

def find_file(pattern, path):
    result = []
    for f in os.listdir(path):
        if fnmatch.fnmatch(f, pattern):
            result.append(f) #os.path.join(root, name))
    return result

def run_cmd(lis,ver=1):
    if ver==1:
        print("********************************")
        print(" ".join(lis))
        print("********************************")
    os.system(" ".join(lis))
    
def one_sample(row):
    sample_name=row[0].split("_")[0]
    if not os.path.isfile(os.path.join(fastas_dir,sample_name+".fasta")):
        r1=os.path.join(fastqs_path,row[0])
        r2=os.path.join(fastqs_path,row[1])
        if os.path.isfile(r1) and os.path.isfile(r2):
            if row[2]!="none":
                long_reads_file=os.path.join(long_reads_path,row[2])
                run_cmd(["unicycler -1",r1,"-2",r2,"-l",long_reads_file,"-o",os.path.join(fastas_dir,sample_name),"-t",str(ncores_per_sample),"--keep 0"])
            else:
                run_cmd(["unicycler -1",r1,"-2",r2,"-o",os.path.join(fastas_dir,sample_name),"-t",str(ncores_per_sample),"--keep 0"])
            
            new_fasta_name=os.path.join(fastas_dir,sample_name+".fasta")
            new_log_name=os.path.join(fastas_dir,sample_name+"_unicycler.log")
            run_cmd(["cp",os.path.join(fastas_dir,sample_name,"assembly.fasta"),new_fasta_name])
            run_cmd(["cp",os.path.join(fastas_dir,sample_name,"unicycler.log"),new_log_name])
            os.system('sed -i "s/>/>'+sample_name+' N/g" '+new_fasta_name)
            os.system('sed -i "s/length/l/g" '+new_fasta_name)
            os.system('sed -i "s/depth/d/g" '+new_fasta_name)
            os.system('sed -i "s/circular=true/cir/g" '+new_fasta_name)
            os.system('sed -i "s/ /_/g" '+new_fasta_name)
            
            run_cmd(["rm -r",os.path.join(fastas_dir,sample_name)])
        else:
            print("One of the fastq files doesn't exist")


###########Command line arguments parshing
fastqs_path=""
long_reads_path=""
R1_pattern=""
out_folder=""
ncores_per_sample="1"
ncores="1"

args=sys.argv
if len(args)>1:
    arguments_file=args[1]

else:
    arguments_file="/home/javiernunezgarcia/AMR_Team_tools/unicycler_hybrid_parallel_arguments_template.args"
    print("Argument file not given or doesn't exist. Please re run with /full/path/to/arguments/file/arguments_file.args")
    #sys.exit()

########Loading other arguments from the arguments file
print('Reading arguments from file: '+arguments_file)
with open(arguments_file,'r') as f:
    for line in f:
        if line[0] not in ['\n',' ','#']:
            print('Loading argument: '+ line.strip())
            exec(line.strip())

############### creating the output directory
if not os.path.exists(out_folder):
    run_cmd(["mkdir -p",out_folder])
run_cmd(['cp ',arguments_file,os.path.join(out_folder,arguments_file.split(os.sep)[-1])])
            
###############checking the samples to run            
print("***** Checking samples to be run")
fastq_R1s=find_file("*"+R1_pattern+"*.fastq.gz", fastqs_path)
R2_pattern=R1_pattern.replace("R1","R2")
fastq_to_process=[]
summary=[]

for fastq_R1 in fastq_R1s:
    fastq_R2=fastq_R1.replace(R1_pattern,R2_pattern)
    sample_name=fastq_R1.split(os.sep)[-1].split("_")[0]
    if os.path.isfile(os.path.join(fastqs_path,fastq_R2)):
        R2_ok="found"
    else:
        R2_ok="not_found"
        print("Missing R2 fastq for file: "+fastq_R1)
        
    if os.path.exists(long_reads_path):    
        long_reads_file=[f for f in listdir(long_reads_path) if sample_name==f[:len(sample_name)] and R2_pattern not in f and R1_pattern not in f]  
        if len(long_reads_file)==1:
            long_reads_file_ok="found"
            long_reads_file=long_reads_file[0]
        elif len(long_reads_file)==0:
            long_reads_file_ok="not_found"
            print("Missing long reads file for file: "+fastq_R1)
        else:
            long_reads_file_ok="several"
            print("More than one long reads file for: "+fastq_R1)
    else:
        print("Only short reads files detected")
        long_reads_file="none"
        long_reads_file_ok="none"
    if R2_ok=="found" and (long_reads_file_ok in ["found","none"] or long_reads_file=="none"):
        fastq_to_process.append([fastq_R1,fastq_R2,long_reads_file])
    
    summary.append([fastq_R1,R2_ok,long_reads_file_ok])
    
print("Check file "+os.path.join(out_folder,"summary.csv"))            

print(fastq_to_process)
writeCSV(os.path.join(out_folder,"summary.csv"),[["R1","R2_status","Long_reads_status"]]+summary)            

     
######################## creating a temporary directory in the home directory with a random name
random.seed(datetime.now())
fastas_dir=os.path.join(str(Path.home()),"uc_fastas"+str(random.randint(0,10000000)))
if not os.path.exists(fastas_dir):
    run_cmd(["mkdir",fastas_dir])

######################## running a pool of samples
pool=mp.Pool(int(ncores))
result=pool.map(one_sample,fastq_to_process)

######################## copying the results into the output directory
run_cmd(["cp -r",fastas_dir,out_folder])

########################deleting the temporary directory
#run_cmd(["rm -r",fastas_dir])




