#!/usr/bin/env python
# coding: utf-8

#AL - the above code is new for the griffin paper version
import argparse
import os
import pandas as pd
import numpy as np

# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Write a script to run copy and reformat "CN" info from ichor')
parser.add_argument('-i','--ichorCNA_output', dest='ichorCNA', 
                   help='results/ichorCNA/ directory',required=True)
args = parser.parse_args()
###########################################################################
ichorCNA=args.ichorCNA.rstrip('/')

try: #make a directory for the ichor output
	os.mkdir('CN_data/')
except FileExistsError:
	pass

for item in sorted(os.listdir(ichorCNA)):
	original_file=ichorCNA+'/'+item+'/'+item+'.correctedDepth.txt'
	if os.path.exists(original_file):
		CN_data=pd.read_csv(original_file, sep='\t')
		CN_data=CN_data[~CN_data['log2_TNratio_corrected'].isnull()]#keep rows that don't have a null log ratio
		output_file='CN_data/'+item+'.correctedDepth_filtered.txt'
		CN_data.to_csv(output_file,sep='\t', index=False)
	else: #for things that are not folders with data
		continue #skip to the next folder
	
		