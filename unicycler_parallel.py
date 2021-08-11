### ASk Javier Nunez for any question
## This script runs unicyler in parallel. It uses the instance hard drive as temporary/auxiliar storage and 
## then it copy the assemblies where they come from in a folder caller "uc_fastas".
## It takes 4 arguments:  
## path to the folder containing the fastq files 
## number of cores to be use for each individual unicycler run
## number of parallel unicycler runs
## pattern to differentiate R1 fastq files. Most common used are "_R1_" or "_R1" or even "R1"
 

import sys,os,os.path
import multiprocessing as mp

def uni_one_sample(r1):
    sample_name=r1.split("_")[0]
    r2=r1.replace(R1_pattern,R1_pattern.replace("R1","R2"))
    r1=os.path.join(fastqs_path,r1)
    r2=os.path.join(fastqs_path,r2)
    print("unicycler -1 "+r1+" -2 "+r2+" -o "+fastas_dir+"/"+sample_name+" -t "+ncores_per_sample+" --keep 1 --vcf")
    os.system("unicycler -1 "+r1+" -2 "+r2+" -o "+fastas_dir+"/"+sample_name+" -t "+ncores_per_sample+" --keep 1 --vcf")
    
    new_fasta_name=fastas_dir+"/"+sample_name+".fasta"
    os.system("cp "+fastas_dir+"/"+sample_name+"/assembly.fasta"+" "+new_fasta_name)
    
    os.system('sed -i "s/>/>'+sample_name+' N/g" '+new_fasta_name)
    os.system('sed -i "s/length/l/g" '+new_fasta_name)
    os.system('sed -i "s/depth/d/g" '+new_fasta_name)
    os.system('sed -i "s/circular=true/cir/g" '+new_fasta_name)
    os.system('sed -i "s/ /_/g" '+new_fasta_name)
    
    os.system("rm -r "+fastas_dir+"/"+sample_name)


###########Command line arguments parshing
fastqs_path=""
R1_pattern=""
out_folder=""
ncores_per_sample="1"
ncores="1"

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

########################################
fastq_R1s=[x for x in os.listdir(fastqs_path) if R1_pattern in x]  

fastas_dir="~/uc_fastas"
os.system("mkdir "+fastas_dir)

pool=mp.Pool(int(ncores))
result=pool.map(uni_one_sample,fastq_R1s)

if not os.path.exists(out_folder):
    os.system("mkdir "+out_folder)
os.system("cp -r "+fastas_dir+" "+out_folder)


