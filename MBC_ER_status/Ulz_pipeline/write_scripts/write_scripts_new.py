#!/usr/bin/env python
# coding: utf-8

#AL - the above code is new for the griffin paper version

import argparse
import os

# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Write scripts to run trim_from_bam_single_end.py, run_tf_analysis_from_bam.py, and scoring pipeline on multiple samples')
parser.add_argument('-i','--input-YAML', dest='in_file', 
                   help='Input YAML file',required=True)
parser.add_argument('-n','--script-name', dest='name',
                   help='name_for_script',required=True)
args = parser.parse_args()
###########################################################################
#get sample list:
with open (args.in_file,'r') as file:
	raw_samples=file.read()

raw_samples=raw_samples.strip('\n\t ').split('\n') #remove empty lines from the end before splitting into lines
raw_samples.remove('samples:')  

samples=[]
for line in raw_samples:
	line=line.strip('\t ')
	if line.startswith('#'):
		continue
	else:
		samples.append(line.strip('\t ').split(': '))
##############################


k = 'trim_align_'+args.name+'.sh' #file name for sh script
current = 'sh_scripts/'+k

try: #make a directory for the sh scripts
	os.mkdir('sh_scripts')
except FileExistsError:
	pass
	
try: #make a directory for the sh scripts
	os.mkdir('trimmed_bams/')
except FileExistsError:
	pass

f = open(current, 'w+')

for sample in samples:
	outpath=sample[1].rsplit('/',1)[1]
	outpath='_60bp.'.join(outpath.rsplit('.',1)) #add 60bp to the bam filename
	
	
	outpath='trimmed_bams/'+outpath
	f.write('time trim_from_bam_single_end.py -i ' + sample[1] + ' -o ' + outpath + ' -l 60 \n')
f.close()

try: #make a directory for the sh scripts
	os.mkdir('results')
except FileExistsError:
	pass

# write run_tf_analysis.py script
run_tf_script='sh_scripts/run_tf_analysis_'+args.name+'.sh'
with open(run_tf_script, 'w+') as f:
	for sample in samples:
		bam_name=sample[1].rsplit('/',1)[1]
		bam_name=bam_name.rsplit('.',1)[0] 

	
		f.write('time run_tf_analyses_from_bam.py -b trimmed_bams/'+bam_name+'_60bp.bam.sorted.bam')
		f.write(' -o '+sample[0])
		f.write(' -norm-file CN_data/'+bam_name+'.correctedDepth_filtered.txt')
		f.write(' -calccov -a tf_gtrd_1000sites\n')

# write run_tf_analysis.py script
run_scoring='sh_scripts/run_scoring_'+args.name+'.sh'
with open(run_scoring, 'w+') as f:
	for sample in samples:
		f.write('time scoring_pipeline.sh '+sample[0])
		f.write(' '+sample[0]+'\n')
		
print ('\nCommands to run scripts:')
print('launchSlurm -p campus-new -b '+current)
print('python write_scripts/reformat_log2_file.py -i <path/to/results/ichorCNA>')
print('launchSlurm -p campus-new -b '+run_tf_script)
print('launchSlurm -p campus-new -b '+run_scoring)
print('\n')

