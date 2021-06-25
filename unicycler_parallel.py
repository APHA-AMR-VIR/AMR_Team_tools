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
    r2=r1.replace(sys.argv[4],sys.argv[4].replace("R1","R2"))
    r1=os.path.join(sys.argv[1],r1)
    r2=os.path.join(sys.argv[1],r2)
    print("unicycler -1 "+r1+" -2 "+r2+" -o "+fastas_dir+"/"+sample_name+" -t "+sys.argv[2]+" --keep 1 --vcf")
    os.system("unicycler -1 "+r1+" -2 "+r2+" -o "+fastas_dir+"/"+sample_name+" -t "+sys.argv[2]+" --keep 1 --vcf")
    
    new_fasta_name=fastas_dir+"/"+sample_name+".fasta"
    os.system("cp "+fastas_dir+"/"+sample_name+"/assembly.fasta"+" "+new_fasta_name)
    
    os.system('sed -i "s/>/>'+sample_name+' N/g" '+new_fasta_name)
    os.system('sed -i "s/length/l/g" '+new_fasta_name)
    os.system('sed -i "s/depth/d/g" '+new_fasta_name)
    os.system('sed -i "s/circular=true/cir/g" '+new_fasta_name)
    os.system('sed -i "s/ /_/g" '+new_fasta_name)
    
    os.system("rm -r "+fastas_dir+"/"+sample_name)

fastq_R1s=[x for x in os.listdir(sys.argv[1]) if sys.argv[4] in x]  

fastas_dir="~/uc_fastas"
os.system("mkdir "+fastas_dir)

pool=mp.Pool(int(sys.argv[3]))
result=pool.map(uni_one_sample,fastq_R1s)

os.system("cp -r "+fastas_dir+" "+sys.argv[1])


