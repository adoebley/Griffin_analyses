#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse
import sys
import os
import pandas as pd
from matplotlib import pyplot as plt


# In[ ]:


print('starting')
sys.stdout.flush()


# In[ ]:


# %load_ext autoreload
# %autoreload 2

# cmd = ['--in_file','../ATAC_nucleosome_profiling/analysis/merged_data/0_0_FC_reformatted.txt',
#        '--name','1_FC',
#        '--status_column','revisions_ER_status_binary',
#       '--iterations','10',
#       '--out_dir','tmp',
#       '--report_interval','2',
#       '--script_path','/fh/fast/ha_g/user/adoebley/projects/griffin_revisions_1/MBC/scripts/']


# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('--in_file', help='/path/to/data_reformatted.txt', required=True)
parser.add_argument('--name', help='name for outputs', required=True)
parser.add_argument('--status_column', help='column containing the status', required=True)
parser.add_argument('--iterations', help='number of iterations', required=True, type = int)
parser.add_argument('--out_dir', help='output directory', required=True)
parser.add_argument('--report_interval', help='how frequently to print progress', default = 50, type = int)
parser.add_argument('--script_path', help='path_to_cancer_detection_logreg_scripts', required = True)


args = parser.parse_args()


# In[ ]:


#import the griffin scripts
sys.path.insert(0, args.script_path)
import MBC_logreg_bootstrap_functions


# In[ ]:


#make output directory
if not os.path.exists(args.out_dir):
    os.mkdir(args.out_dir)


# In[ ]:


print('importing data')
data,features  = MBC_logreg_bootstrap_functions.import_data(args.in_file,args.status_column,0.05, 0.1)
sys.stdout.flush()

print('all samples')
print(data[args.status_column].value_counts(dropna = False))

print('first samples')
print(data[data['first_passing_sample']==1][args.status_column].value_counts(dropna = False))
sys.stdout.flush()

print('running '+str(args.iterations)+' logreg bootstrap iterations')
probabilities,c_vals,coefs = MBC_logreg_bootstrap_functions.run_bootstrap(data,args.iterations,features,args.report_interval)    
sys.stdout.flush()

#export results

    
probabilities.to_csv(args.out_dir+'/'+args.name+'.probabilities.txt', sep='\t', float_format='%.5f')
pd.Series(c_vals).to_csv(args.out_dir+'/'+args.name+'.c_values.txt', sep='\t', header = False, index=False)
coefs.to_csv(args.out_dir+'/'+args.name+'.coefs.txt', sep='\t', float_format='%.5f')

plt.hist([str(m) for m in sorted(c_vals)])
plt.savefig(args.out_dir+'/'+args.name+'.cvals.pdf')   

print('Getting accuracy and AUC')
sys.stdout.flush()
accuracies,accuracy_CIs,AUCs,AUC_CIs = MBC_logreg_bootstrap_functions.get_accuracy_AUC(probabilities,data,args.iterations)

accuracies.to_csv(args.out_dir+'/'+args.name+'.accuracy.txt', sep='\t', index = False, float_format='%.5f')
accuracy_CIs.to_csv(args.out_dir+'/'+args.name+'.accuracy_CI.txt', sep='\t', float_format = '%.5f')
AUCs.to_csv(args.out_dir+'/'+args.name+'.AUC.txt', sep='\t', index = False, float_format='%.5f')
AUC_CIs.to_csv(args.out_dir+'/'+args.name+'.AUC_CI.txt', sep='\t', float_format = '%.5f')

print('Done')


# In[ ]:




