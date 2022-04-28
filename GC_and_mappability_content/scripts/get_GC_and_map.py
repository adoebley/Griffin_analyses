#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#import stuff 
import pandas as pd
import numpy as np

from scipy.signal import savgol_filter

import os
import sys
import time

import pysam
import pybedtools
import pyBigWig

import argparse


# In[ ]:


# from matplotlib import pyplot as plt
# %matplotlib inline

# #define paths and parameters
# TF_name = 'LYL1'
# site_list='../sites/TFBS/10000_unfiltered_sites_CIS_BP_v2/'+TF_name+'.hg38.10000.txt'

# tmp_dir = 'tmp'
# results_dir = 'results'


# ref_seq_path='/fh/fast/ha_g/grp/reference/GRCh38/GRCh38.fa'
# mappability_bw='../genome/k100.Umap.MultiTrackMappability.hg38.bw'
# chrom_sizes_path = '/fh/fast/ha_g/grp/reference/GRCh38/hg38.standard.chrom.sizes'

# #define the columns in the TFBS files
# chrom_col='Chrom'
# pos_col='position'
# chroms = ['chr'+str(m) for m in np.arange(1,23)]


# window_start=-1000
# window_end=1000

# encode_exclude = '../genome/encode_unified_GRCh38_exclusion_list.bed'
# centromere_path = '../genome/hg38_centromeres.bed'
# gap_path = '../genome/hg38_gaps.bed'
# patch_path = '../genome/hg38_fix_patches.bed'
# alternative_haplotype_path = '../genome/hg38_alternative_haplotypes.bed'

# exclude_paths = [encode_exclude,centromere_path,gap_path,patch_path,alternative_haplotype_path]
# del(encode_exclude,centromere_path,gap_path,patch_path)


# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('--TF_name', help='name of TF', required=True)
parser.add_argument('--site_list', help='list of sites', required=True)

parser.add_argument('--tmp_dir', help = 'directory for temporary outputs', required=True)
parser.add_argument('--results_dir', help = 'directory for results', required=True)

parser.add_argument('--reference_genome',help = 'path to the reference genome',required=True)
parser.add_argument('--mappability_bw',help = 'bigWig file of genome wide mappability scores',required=True)
parser.add_argument('--chrom_sizes_path', help='path to chrom sizes file', required=True)

parser.add_argument('--chrom_column',help='name of column containing chromosome number', required = True)
parser.add_argument('--position_column',help='name of column containing chromosome position', required = True)
parser.add_argument('--chroms', help='chroms to include when selecting sites', nargs='*', required = True)

parser.add_argument('--save_window',help='start and end of the window around each site',nargs=2, type=int, required = True)

parser.add_argument('--exclude_paths', help='path to bed files of regions to filter out (excluded regions, centromeres, gaps, patches, alternative haplotypes), or "none" to not exclude any regions', required=True, nargs = '*')


args = parser.parse_args()


TF_name = args.TF_name
site_list = args.site_list

tmp_dir = args.tmp_dir 
results_dir = args.results_dir 

ref_seq_path = args.reference_genome
mappability_bw = args.mappability_bw
chrom_sizes_path = args.chrom_sizes_path

chrom_col = args.chrom_column
pos_col=args.position_column
chroms = args.chroms

window_start,window_end = args.save_window

exclude_paths = args.exclude_paths


# In[ ]:


#open ref seq
ref_seq=pysam.FastaFile(ref_seq_path)

#get the window
window_len=window_end-window_start
window_values=np.arange(window_start,window_end)

print('window_length:',window_len)
sys.stdout.flush()


# In[ ]:


#snakemake should create these folders, but if not using the snakemake, this is needed
if not os.path.exists(tmp_dir): 
    os.mkdir(tmp_dir)

tmp_pybedtools = tmp_dir
pybedtools.set_tempdir(tmp_pybedtools)

tmp_bigWig = tmp_dir


# In[ ]:


merged_exclude_regions = pybedtools.BedTool('\n', from_string=True)

#create an empty bed file
excluded_regions_bw = pyBigWig.open(tmp_bigWig+"/excluded_regions.bw", "w")
chrom_sizes = pd.read_csv(chrom_sizes_path, sep='\t', header=None)
chrom_sizes = chrom_sizes[chrom_sizes[0].isin(chroms)]
excluded_regions_bw.addHeader([(a,b) for a,b in chrom_sizes.values])

for path in exclude_paths:
    print('excluding:',path)
    current_regions = pybedtools.BedTool(path)
    merged_exclude_regions = merged_exclude_regions.cat(current_regions)    
    del(current_regions)
merged_exclude_regions = merged_exclude_regions.to_dataframe()
merged_exclude_regions = merged_exclude_regions[merged_exclude_regions['chrom'].isin(chroms)]
pybedtools.cleanup()
excluded_regions_bw.addEntries(list(merged_exclude_regions['chrom']), list(merged_exclude_regions['start']), ends = list(merged_exclude_regions['end']), values = [1.0 for m in range(len(merged_exclude_regions))])  

excluded_regions_bw.close()
sys.stdout.flush()


# In[ ]:


sites=pd.read_csv(site_list, sep='\t')

print(len(sites))
print(sites.head())


# In[ ]:


#fetch zero mappability
print('Fetching mappability')
start_time = time.time()

#make a np array for holding data
mappability_mask=pd.DataFrame(columns = window_values)
bw = pyBigWig.open(mappability_bw, "r")

for i in range(len(sites)):
    if i%1000==0:
        print(i,time.time()-start_time)
        sys.stdout.flush()
    #get the location of the site (chromosome and center of site)
    chrom = sites.iloc[i][chrom_col]
    position = sites.iloc[i][pos_col]
    
    
    fetched=bw.values(chrom,position+window_start,position+window_end, numpy = True)

    fetched = pd.Series(fetched, index = window_values)
    mappability_mask = mappability_mask.append(fetched, ignore_index=True)
    
bw.close()


# In[ ]:


#fetch excluded sites
print('Fetching excluded sites')

#make a np array for holding data
exclusion_mask=pd.DataFrame(columns = window_values)
bw = pyBigWig.open(tmp_bigWig+"/excluded_regions.bw", "r")

for i in range(len(sites)):
    if i%1000==0:
        print(i,time.time()-start_time)
        sys.stdout.flush()
    #get the location of the site (chromosome and center of site)
    chrom = sites.iloc[i][chrom_col]
    position = sites.iloc[i][pos_col]
    
    fetched=bw.values(chrom,position+window_start,position+window_end, numpy = True)

    fetched = pd.Series(fetched, index = window_values)
    exclusion_mask = exclusion_mask.append(fetched, ignore_index=True)
    
bw.close()


# In[ ]:


print('Fetching GC content')
#make a np array for holding data
site_values=pd.DataFrame(columns = window_values)

for i in range(len(sites)):
    if i%1000==0:
        print(i,time.time()-start_time)
        sys.stdout.flush()
    #get the location of the site (chromosome and center of site)
    chrom = sites.iloc[i][chrom_col]
    position = sites.iloc[i][pos_col]
    
    fetched=ref_seq.fetch(reference=chrom,start=position+window_start,end=position+window_end)
    fetched = np.array(list(fetched.upper()))
    fetched[np.isin(fetched, ['A','T','W'])] = 0
    fetched[np.isin(fetched, ['C','G','S'])] = 1
    rng = np.random.default_rng(position)
    fetched[np.isin(fetched, ['N','R','Y','K','M','B','D','H','V'])] = rng.integers(2, size=len(fetched[np.isin(fetched, ['N','R','Y','K','M','B','D','H','V'])])) #random integer in range(2) (i.e. 0 or 1)
    fetched = fetched.astype(float)
    fetched = pd.Series(fetched, index = window_values)
    site_values = site_values.append(fetched, ignore_index=True)
    
    


# In[ ]:


#exclude zero mappability
results = np.where(mappability_mask>0,site_values,np.nan)#only retain positions where the mappability is >0
results = pd.DataFrame(results)
results.columns = window_values

#exclude excluded regions
results = np.where(exclusion_mask>0,np.nan,results)
results = pd.DataFrame(results)
results.columns = window_values

#savgol smoothing for GC values 
all_smoothed_sites = savgol_filter(results,165,0)
all_smoothed_sites = pd.DataFrame(all_smoothed_sites)
all_smoothed_sites.columns = window_values

#add ID column
all_smoothed_sites['unique_ID'] = sites[chrom_col]+'_'+sites[pos_col].astype(str)
all_smoothed_sites = all_smoothed_sites[['unique_ID']+list(window_values)]
all_smoothed_sites = all_smoothed_sites.set_index('unique_ID')


# In[ ]:


#exclude regions from mappability values
mappability_results = np.where(exclusion_mask>0,np.nan,mappability_mask)
mappability_results = pd.DataFrame(mappability_results)
mappability_results.columns = window_values

#these don't get smoothed 

#add ID column
mappability_results['unique_ID'] = sites[chrom_col]+'_'+sites[pos_col].astype(str)
mappability_results = mappability_results[['unique_ID']+list(window_values)]
mappability_results = mappability_results.set_index('unique_ID')


# In[ ]:


all_smoothed_sites.to_csv(results_dir+'/'+TF_name+'.smoothed_GC_content.tsv', sep='\t', float_format='%.2f')
mappability_results.to_csv(results_dir+'/'+TF_name+'.mappability.tsv', sep='\t')


# In[ ]:


pybedtools.cleanup('all')
os.remove(tmp_bigWig+"/excluded_regions.bw")


# In[ ]:





# In[ ]:


# #GC content plot
# plt.plot(window_values,np.nanmean(all_smoothed_sites, axis = 0), label = 'smoothed GC content') 
# plt.fill_between(window_values, np.nanpercentile(all_smoothed_sites,25, axis = 0), np.nanpercentile(all_smoothed_sites,75, axis = 0), alpha = 0.25, label = 'interquartile range')

# plt.legend()
# plt.xlabel('distance from TFBS')


# In[ ]:


# #mappability plot
# plt.plot(window_values,np.nanmean(mappability_results, axis = 0), label = 'mean mappability') 
# # plt.fill_between(window_values, np.nanpercentile(mappability_results,25, axis = 0), np.nanpercentile(all_smoothed_sites,75, axis = 0), alpha = 0.25, label = 'interquartile range')
# plt.ylim(.9895,1.0005)

# plt.legend()
# plt.xlabel('distance from TFBS')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




