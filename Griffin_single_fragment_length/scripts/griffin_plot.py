#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import time
import sys
import yaml
import argparse


# In[ ]:


# #for ipynb
# %matplotlib inline

# in_dir = '../snakemake/griffin_nucleosome_profiling/results/'
# samples_yaml = '../snakemake/griffin_nucleosome_profiling/config/samples.GC.yaml'

# save_window = [-500,500]
# step = 15

# individual = 'True'
# out_dir = 'tmp'


# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('--in_dir', help='path/to/results/', required=True)
parser.add_argument('--samples_yaml', help='samples.GC.yaml', required=True)

parser.add_argument('--save_window',help='start and end of window to be plotted',nargs=2, type=int, default=(-1000,1000))
parser.add_argument('--step',help='step size when calculating coverage', type=int, default=5)

parser.add_argument('--individual',help='if individual sites were saved in previous steps. (True/False)',default='False')
parser.add_argument('--out_dir',help='folder for results',required=True)

args = parser.parse_args()

in_dir = args.in_dir
samples_yaml = args.samples_yaml

save_window=args.save_window
step = args.step

individual = args.individual
out_dir = args.out_dir


# In[ ]:


save_window=[int(np.ceil(save_window[0]/step)*step),int(np.floor(save_window[1]/step)*step)] #round to the nearest step inside the window
save_columns = np.arange(save_window[0],save_window[1],step)
str_save_columns = [str(m) for m in save_columns]
print('save_window rounded to step:',save_window)


# In[ ]:


with open(samples_yaml,'r') as f:
    samples = yaml.safe_load(f)
    
samples = samples['samples']
samples = list(samples.keys())


# In[ ]:


#dict to hold results grouped by correction type
print('Importing data')
results_dict = {'uncorrected': pd.DataFrame(),
                'GC_corrected': pd.DataFrame(),
                'GC_map_corrected': pd.DataFrame()}
#import
for sample in samples:
    print(sample)
    for key in results_dict.keys():
        current_file = in_dir+'/'+sample+'/'+sample+'.'+key+'.coverage.tsv'
        current = pd.read_csv(current_file, sep='\t')
        if individual.lower()=='true':
            current = current.groupby('site_name')[str_save_columns].mean()
            current = current.reset_index() #retain site_name
            current['sample'] = sample
        results_dict[key] = results_dict[key].append(current, ignore_index=True)

site_names = results_dict['uncorrected']['site_name'].unique()


# In[ ]:


#if info about individual sites was kept, the averaging process can take quite a while. Save for later use. 
if individual.lower()=='true':
    for i,key in enumerate(results_dict.keys()):
        data = results_dict[key].copy()

        data.to_csv(out_dir+'/plots/'+key+'.mean_data.txt', sep='\t', index=False)
    


# In[ ]:


#generate plots
for j,site_name in enumerate(site_names):
    fig,axes = plt.subplots(1,3,figsize=(12,3.5), sharey = 'row')
    for i,key in enumerate(results_dict.keys()):
        data = results_dict[key].copy()
        ax = axes[i]
        for sample in data['sample'].unique():
            current = data[(data['sample']==sample) & (data['site_name']==site_name)]
            ax.plot(save_columns, current[str_save_columns].T, label=sample)
            ax.tick_params(labelleft=True)
        ax.set_title(site_name+' '+key)
        ax.set_xlabel('distance from site')
    
    axes[0].set_ylabel('normalized coverage')
    
    if len(data['sample'].unique())<15:
        axes[2].legend(bbox_to_anchor=[1,1],loc = 'upper left')
    else:
        axes[2].legend(bbox_to_anchor=[1,1],loc = 'upper left',ncol=2)

    fig.tight_layout()
    plt.savefig(out_dir+'/plots/'+site_name+'.pdf')
    plt.close('all')
    if j%20==0:
        print(j,site_name)
        sys.stdout.flush()


# In[ ]:





# In[ ]:




