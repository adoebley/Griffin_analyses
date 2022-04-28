#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pysam
import os
#import pybedtools #not used
import pandas as pd
import numpy as np
import time
import argparse
import sys
from matplotlib import pyplot as plt


# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('--bam_file_name', help='sample name (does not need to match actual file name)', required=True)

parser.add_argument('--mappable_name', help='name of mappable regions file (with .bed removed)', required=True)

parser.add_argument('--genome_GC_frequency',help='folder containing GC counts in the reference sequence (made by generate_reference_files.snakemake)',required=True)

parser.add_argument('--out_dir',help='folder for GC bias results',required=True)

parser.add_argument('--size_range',help='range of read sizes to be included',nargs=2, type=int, required=True)
parser.add_argument('--fragment_length',help='most common fragment size',type=int, required=True)

args = parser.parse_args()

bam_file_name = args.bam_file_name
mappable_name=args.mappable_name
genome_GC_frequency = args.genome_GC_frequency
out_dir = args.out_dir
size_range = args.size_range
fragment_length = args.fragment_length


# In[ ]:


print('arguments provided:')

print('\tbam_file_name = "'+bam_file_name+'"')
print('\tmappable_name = "'+mappable_name+'"')

print('\tgenome_GC_frequency = "'+genome_GC_frequency+'"')
out_dir = out_dir.rstrip('/')
print('\tout_dir = "'+out_dir+'"')

print('\tsize_range = '+str(size_range))
print('\tfragment_length = '+str(fragment_length))


# In[ ]:


#For now I'm going to keep the smoothing bin size as a set variable
GC_smoothing_step = 20


# In[ ]:


#input is the out_file from the previous step
in_file = out_dir +'/GC_counts/'+ bam_file_name+'.GC_counts.txt'
print('in_file:',in_file)

#output is smoothed version
smoothed_out_file = out_dir +'/GC_bias/'+ bam_file_name+'.GC_bias.txt'

#plot files
plot_file1 = out_dir +'/GC_plots/'+ bam_file_name+'.GC_bias.pdf'

print('out_file:',smoothed_out_file)
sys.stdout.flush()


# In[ ]:


#create output folders if needed
if not os.path.exists(out_dir +'/GC_plots/'):
    os.mkdir(out_dir +'/GC_plots/')
if not os.path.exists(out_dir +'/GC_bias/'):
    os.mkdir(out_dir +'/GC_bias/')


# In[ ]:


#import the GC info from the genome
current_path = genome_GC_frequency+'/'+mappable_name+'.'+str(fragment_length)+'bp.GC_frequency.txt'
GC_freq = pd.read_csv(current_path,sep='\t')
    
GC_freq['GC_content']=GC_freq['num_GC']/GC_freq['length']
GC_freq = GC_freq.sort_values(by=['GC_content','length']).reset_index(drop=True)


# In[ ]:


#import GC counts from the sample
GC_df = pd.read_csv(in_file, sep='\t')

GC_df['GC_content']=GC_df['num_GC']/GC_df['length']
GC_df = GC_df.sort_values(by=['GC_content','length']).reset_index(drop=True)


# In[ ]:


#calculate the GC_bias

new_df = GC_df[GC_df['length']==fragment_length].copy().reset_index(drop=True)
current_freq = GC_freq[GC_freq['length']==fragment_length].copy().reset_index(drop=True)

#save the frequency of each GC content in the genome
new_df['number_of_positions']=current_freq['number_of_fragments']

#calculate the GC bias
current_bias = new_df['number_of_fragments']/new_df['number_of_positions']    
new_df['GC_bias'] = current_bias

#normalize to a mean of 1 for each fragment length(compute GC bias does this same thing)
new_df['GC_bias'] = new_df['GC_bias']/np.nanmean(new_df['GC_bias'])

#print(length,len(new_df['GC_bias']),np.nanmean(new_df['GC_bias']))
    
new_df = new_df.sort_values(by=['GC_content','length']).reset_index(drop=True)


# In[ ]:


def median_smoothing(current,fraction):
    bin_size=int(len(current)*fraction)
    #if bin_size<10:
        #bin_size=10
    medians = []

    for i in range(len(current)):
        start = int(i-bin_size/2)
        end = int(i+bin_size/2)
        #if the bin starts before the beginning, just take the first bin
        if start<0:
            start=0
            end=bin_size
        #if the bin extends beyond the end, take the last bin
        if end>=len(current):
            start=len(current)-bin_size
            end=len(current)
        current_median = np.nanmedian(current['GC_bias'].iloc[start:end])
        medians.append(current_median)
    return(medians)


# In[ ]:



#perform smoothing
fit = median_smoothing(new_df,.05)  
new_df['smoothed_GC_bias']=fit
    
    
#get rid of values for GC contents that are never observed
new_df['smoothed_GC_bias'] = np.where(new_df['number_of_positions']==0,np.nan,new_df['smoothed_GC_bias'])
    
#normalize to a mean of 1
new_df['smoothed_GC_bias'] = new_df['smoothed_GC_bias']/np.nanmean(new_df['smoothed_GC_bias'])
    


# In[ ]:


#export results
new_df.to_csv(smoothed_out_file,sep='\t',index=False)


# In[ ]:


fig, axes = plt.subplots(1,2, figsize = (10,4))

axes[0].plot(new_df['GC_content'], new_df['number_of_positions']/np.nanmean(new_df['number_of_positions']),label = 'number_of_positions')
axes[0].plot(new_df['GC_content'], new_df['number_of_fragments']/np.nanmean(new_df['number_of_fragments']),label = 'number_of_fragments')
axes[0].set_ylabel('Frequency')
axes[0].set_xlabel('GC_content')

axes[0].legend(bbox_to_anchor = [1,1], loc = 'upper left')

axes[1].plot(new_df['GC_content'], new_df['GC_bias'],label = 'GC_bias')
axes[1].plot(new_df['GC_content'], new_df['smoothed_GC_bias'],label = 'smoothed_GC_bias')
axes[1].set_ylabel('GC_bias')
axes[1].set_xlabel('GC_content')
axes[1].legend(bbox_to_anchor = [1,1], loc = 'upper left')

fig.suptitle(bam_file_name+'\n'+mappable_name + '\n' + str(fragment_length)+'bp')
fig.tight_layout()
fig.savefig(plot_file1)
plt.close()


# In[ ]:




