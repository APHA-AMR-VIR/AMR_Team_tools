### ASk Javier Nunez for any question
## This script runs prokka in parallel. It uses the instance hard drive as temporary/auxiliar storage and 
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
    
def one_sample(fasta_name):
    fasta=os.path.join(fastas_path,fasta_name)
    sample_name=fasta_name.split(".")[0].split("_")[0]
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
    print(fasta)
    print(sample_name)
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
    run_cmd(['prokka','--outdir',os.path.join(aux_dir,sample_name),'--locustag',sample_name,'--prefix',sample_name+'.pk',fasta,'--cpus',str(ncores_per_sample),'--force'])
    #if not os.path.isfile(os.path.join(aux_dir,sample_name+'.pk')):
    #    return("Prokka didn't work for fasta file "+fasta)
    #else:
    #    return("Prokka done for fasta file "+fasta)



#########Command line arguments parshing
fastas_path=""
extension_fastas=""
out_folder=""
ncores_per_sample="1"
ncores="1"

args=sys.argv
if len(args)>1:
    arguments_file=args[1]

else:
    arguments_file="/home/javiernunezgarcia/AMR_Team_tools/AMR_prokka_arguments_template.args"
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

fastas=[fil for fil in os.listdir(fastas_path) if fil.split(".")[-1] in extension_fastas]
print(str(len(fastas))+" fasta files to be processed")
value=input('Happy to go (y/n)?\n')
if value!='y':
    sys.exit() 

####################### creating a temporary directory in the home directory with a random name
random.seed(datetime.now())
aux_dir=os.path.join(str(Path.home()),"prokka_"+str(random.randint(0,10000000)))
if not os.path.exists(aux_dir):
    run_cmd(["mkdir",aux_dir])

######################## running a pool of samples
#for fasta in fastas:
#    one_sample(fasta)
pool=mp.Pool(int(ncores))
result=pool.map(one_sample,fastas)

######################## copying the results into the output directory
run_cmd(["cp -r",aux_dir+'/*',out_folder+'/.'])

########################deleting the temporary directory
run_cmd(["rm -r",aux_dir])




