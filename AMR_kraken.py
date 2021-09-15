### ASk Javier Nunez for any question
## This script runs kraken. It uses the instance hard drive as temporary/auxiliar storage called kraken_xxxxxxx
## which is deleted after the results are copied into the path of your choice
## Check the arguments file associted (same name) for info about the arguments you need to pass
 

import sys,os,os.path,fnmatch,csv
import random
from pathlib import Path
from datetime import datetime
import psutil



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
    r1=os.path.join(fastqs_path,row[0])
    r2=os.path.join(fastqs_path,row[1])
    db=os.path.join(db_dir,database_folder.split(os.sep)[-1])
    report_file=os.path.join(aux_dir,sample_name+'_kraken_report')
    if os.path.isfile(r1) and os.path.isfile(r2):
        run_cmd(['kraken2 --db',db,'--gzip-compressed','--paired',r1,r2,'--threads',str(threads),'--use-names', \
                 '--report',report_file,'--report-zero-counts','--output -']) #os.path.join(aux_dir,sample_name+'_kraken')])
        return("run succesfully? "+str(os.path.isfile(report_file)))
    else:
        print("One of the fastq files doesn't exist")
        return("Wrong_fasta")
    


#concatenating the kraken reports from multiple file using a dictionary
def kraken_cat_report(patho):
    filelist=[os.path.join(patho,x) for x in os.listdir(patho) if "_kraken_report" in x]
    tab=[["Taxa","ID","Highest Species","Highest reads","sum of Reads","% in top species"]]
    taxa_dict=dict()
    for file in filelist:
        sample_name=file.split(os.sep)[-1].split("_")[0]
        highest_val=0
        highest_taxa=""
        sum_reads=0
        taxa_dict[sample_name]=dict()
        openfi=open(file, 'r')
        for line in openfi:
            fields=line.split('\t')
            if (fields[3]=="S"):
                taxa=fields[5].strip()
                reads=int(fields[2])
                taxa_dict[sample_name][taxa]=reads
                if reads>highest_val:
                    highest_val=reads
                    highest_taxa=taxa
                sum_reads=sum_reads+reads
        tab.append([file,sample_name,highest_taxa,highest_val,sum_reads,round(100*highest_val/sum_reads,2)])
    
    samples=[x[1] for x in tab[1:]]
    taxas=[taxa_dict[sample].keys() for sample in samples]
    assert (len(set([len(t) for t in taxas]))==1), "No combination of reports done. Kraken reports differing in the number of lines"
    taxas=list(taxa_dict[samples[0]].keys())
    tab[0]=tab[0]+taxas
    for i in range(1,len(tab)):
        for taxa in taxas:
            tab[i].append(taxa_dict[tab[i][1]][taxa])
    tab=list(map(list, zip(*tab)))
    
    tab_filter=tab[:6]
    for row in tab[6:]:
        if sum(row[1:])>0:
           tab_filter.append(row)
   
    writeCSV(os.path.join(patho,"combined_kraken_reports.csv"),tab_filter)
    return([len(samples),len(taxas)])
        
####################################################
###########Command line arguments parshing
database_folder=""
fastqs_path=""
R1_pattern=""
out_folder=""
ncores="1"
threads="1"

args=sys.argv
if len(args)>1:
    arguments_file=args[1]

else:
    arguments_file="/home/javiernunezgarcia/mnt/fsx-016/javi_temp/AMR_Team_tools/kraken_parallel_arguments_template.ar"
    print("Argument file not given or doesn't exist. Please re run with /full/path/to/arguments/file/arguments_file.args")
    #sys.exit()

########Loading other arguments from the arguments file
print('Reading arguments from file: '+arguments_file)
with open(arguments_file,'r') as f:
    for line in f:
        if line[0] not in ['\n',' ','#']:
            print('Loading argument: '+ line.strip())
            exec(line.strip())

############### loading database in memory
db_dir=os.path.join(str(Path.home()),'mnt','kraken_db')
mem=psutil.virtual_memory().available/(1024*1024*1024) 
value=input('Is the database already loaded into the ramdisk drive '+db_dir+'  (y/n)?\n')
if value=='y':
    assert (mem>55), "Not enough memory to run. You need at least 55GB, only "+str(mem)+" available" 
else:
    assert (mem>115), "Not enough memory to run. You need at least 115GB, only "+str(mem)+" available"
    run_cmd(['mkdir',db_dir])
    run_cmd(['sudo mount -t tmpfs -o size=55000m tmpfs',db_dir])
    print("Loading database into randisk, it could take up to a couple of hours")
    run_cmd(['cp -r',database_folder,db_dir])



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
    if R2_ok=="found":
        fastq_to_process.append([fastq_R1,fastq_R2])
    
    summary.append([fastq_R1,R2_ok])
    
print("Check file "+os.path.join(out_folder,"pre_summary_findings.csv"))            

print(fastq_to_process)
writeCSV(os.path.join(out_folder,"pre_summary_findings.csv"),[["R1","R2_status"]]+summary)            

     
######################## creating a temporary directory in the home directory with a random name
random.seed(datetime.now())
aux_dir=os.path.join(str(Path.home()),"kraken_"+str(random.randint(0,10000000)))
if not os.path.exists(aux_dir):
    run_cmd(["mkdir",aux_dir])

######################## running the samples sequentially
for i in range(len(fastq_to_process)):
    fastq_to_process[i].append(one_sample(fastq_to_process[i]))

writeCSV(os.path.join(aux_dir,"post_results_check_list.csv"),[["R1","R2","finish_ok"]]+fastq_to_process) 

######################### concatenating reports
[n_samples,n_taxas]=kraken_cat_report(aux_dir)
print(str(n_samples)+" concatenated for "+str(n_taxas)+" taxas")
######################## copying the results into the output directory

run_cmd(["rm",aux_dir+'/*_kraken_report']) 

run_cmd(["cp -r",aux_dir+'/*',out_folder+'/.'])

########################deleting the temporary directory
run_cmd(["rm -r",aux_dir])


