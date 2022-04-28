#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse
import sys
import os
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np


# In[ ]:


print('starting')
sys.stdout.flush()


# In[ ]:


# %matplotlib inline

# cmd = ['--in_file','../number_of_sites_analysis/merged_data/30000-sites_reformatted.txt',
#        '--name','30000',
#       '--iterations','1000',
#       '--out_dir','tmp/',
#       '--report_interval','50',
#        '--fraction_variance','.8',
#       '--tfx_groups','0.03','0.05']


# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('--in_file', help='/path/to/data_reformatted.txt', required=True)
parser.add_argument('--name', help='name for outputs', required=True)
parser.add_argument('--iterations', help='number of iterations', required=True, type = int)
parser.add_argument('--out_dir', help='output directory', required=True)
parser.add_argument('--report_interval', help='how frequently to print progress', default = 50, type = int)
parser.add_argument('--fraction_variance', help='what fraction of the variance should be explained by the PCs used in the model', type = float)
parser.add_argument('--tfx_groups', help='max value for each tfx group, in ascending order', type = float, nargs = '*', required = True)

args = parser.parse_args()


# In[ ]:


#create directory to export results
if not os.path.exists(args.out_dir):
    os.mkdir(args.out_dir)


# In[ ]:


def import_data(in_file):
    import pandas as pd
    import numpy as np
    from sklearn.preprocessing import StandardScaler
    
    data = pd.read_csv(in_file, sep='\t')
    data = data.set_index('sample')

    #get features and exclude all other columns
    features = data.columns[(data.columns.str.startswith('central_cov')) | (data.columns.str.startswith('mean_cov')) | (data.columns.str.startswith('amplitude')) | (data.columns.str.startswith('Ulz'))]
    print('Features',len(features))

    data = data.sort_index()

    print('Total samples:',len(data))

    #scale data
    scaler = StandardScaler()
    scaler.fit(data[features])
    data[features] = scaler.transform(data[features])
    data[features].mean()
    
    #add tumor fraction groups
    data['tfx_group'] = 'none'
    a = 0
    for b in args.tfx_groups:
        tfx_group_name = str(a)+'-'+str(b)+'TFx'
        data['tfx_group'] = np.where((data['status']==1) & (data['tumor_fraction']>=a) & (data['tumor_fraction']<b),tfx_group_name,data['tfx_group'])
        a=b
    #if group maxes don't go all the way to 1, add a group > max val
    if b<1:
        tfx_group_name = '>'+str(b)+'TFx'
        data['tfx_group'] = np.where((data['status']==1) & (data['status']==1) & (data['tumor_fraction']>=b),tfx_group_name,data['tfx_group'])
    #specify tfx group for healthy donors
    data['tfx_group'] = np.where((data['status']==0),'Healthy',data['tfx_group'])
    
    return(data,features)


# In[ ]:


#bootstrapping
def run_bootstrap_with_PCA(data,iterations,features,report_interval,hyperparameters):
    import time
    import sys
    import pandas as pd
    import numpy as np
    from sklearn.model_selection import StratifiedKFold
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GridSearchCV
    from matplotlib import pyplot as plt
    from sklearn.decomposition import PCA

    start_time = time.time()

    probabilities = pd.DataFrame(index=data.index)
    c_vals = []
    coefs = pd.DataFrame(index=features)
    num_pcs = []
    train_indexes = []
    
    # Loop for each iteration
    for i in range(iterations):
            
        #bootstrap a training set with replacement
        X_train = data.sample(len(data), replace = True, random_state = i+100)[features]
        y_train = data.sample(len(data), replace = True, random_state = i+100)['status']
        
        #the test set is all samples that aren't seen in the training data
        X_test = data[~(data.index.isin(X_train.index))][features]
        y_test = data[~(data.index.isin(X_train.index))]['status']
        
        #print(len(X_train),len(X_train.index.unique()),len(X_test))
        
        #perform PCA on the training set
        n_components = min(len(features), len(X_train))
        pca = PCA(n_components=n_components, svd_solver='randomized', random_state = 100)
        PCs = pca.fit_transform(X_train[features])
        principal_components = pd.DataFrame(data = PCs, columns = ['PC_'+str(m) for m in np.arange(n_components)], index = X_train.index)
        
        #find the principle components that make up 80% of the varience
        for j in range(len(pca.explained_variance_ratio_)):
            current_sum = pca.explained_variance_ratio_[:j].sum()
            if current_sum>=args.fraction_variance:
                break
        #print('number of components:',j)
        pca_features = ['PC_'+str(m) for m in np.arange(0,j)]
        
        #apply to the test data
        test_PCs = pca.transform(X_test[features])
        test_principal_components = pd.DataFrame(data = test_PCs , columns = ['PC_'+str(m) for m in np.arange(n_components)], index = X_test.index)
        
        X_train = principal_components[pca_features]
        X_test = test_principal_components[pca_features]
        
        #10 fold cross validation on the training set
        cv = StratifiedKFold(n_splits=10, shuffle=True, random_state = i+100) 

        model = LogisticRegression(class_weight='balanced', max_iter=500, solver = 'liblinear')
        search = GridSearchCV(estimator=model, param_grid=hyperparameters, cv=cv, n_jobs = 1)
        search.fit(X_train, y_train)
        best_C = search.best_params_['C']

        ##train a new model on the full training dataset (is this the same as refit...?)
        model = LogisticRegression(class_weight='balanced', max_iter=500, C=best_C, solver = 'liblinear')
        model.fit(X_train, y_train)

        #predict the test data
        pred = model.predict(X_test)
        prob = model.predict_proba(X_test)

        #save results
        probabilities[i] = pd.Series(prob[:,1], index = X_test.index)
        c_vals.append(best_C)
        coefs[i] = pd.Series(model.coef_[0], index = pca_features)
        num_pcs.append(j)
     
        train_indexes.append(list(X_train.index))
        
        if i%report_interval==0:
            print('iteration:',i, ', time (sec):',np.round(time.time()-start_time,2),'num_pcs:',j)
        if i%20==0:
            #prevent dfs from becoming too fragmented
            probabilities = probabilities.copy()
            coefs = coefs.copy()   
            sys.stdout.flush()

    probabilities = probabilities.merge(data[['status']], left_index=True, right_index=True)

    return(probabilities,c_vals,coefs,num_pcs,train_indexes)


# In[ ]:


def get_AUC(probabilities,data,iterations):
    #get AUC and accuracy for each bootstrap
    from sklearn.metrics import roc_curve,auc
    import pandas as pd
    import numpy as np

    AUCs = pd.DataFrame()

    probabilities = probabilities.merge(data[['tumor_fraction','sample_type','Stage','tfx_group']], left_index=True, right_index=True)
    
    for i in range(iterations):
        current_dict = {}
        current = probabilities[~(probabilities[i].isnull())][['status','tumor_fraction','sample_type','Stage','tfx_group',i]].copy()

        #overall accuracy and AUC
        group = 'overall'
        fpr,tpr,_ = roc_curve(current['status'],current[i])
        AUC = auc(fpr,tpr)
        current_dict[group] = AUC
        del(AUC,group,fpr,tpr)

        #separate out the healthy samples to be used in every AUC
        healthy_df = current[current['status']==0]
        cancer_df = current[current['status']==1]
        del(current)
        
        for group,df in cancer_df.groupby('sample_type'):
            if group == 'Duodenal_Cancer':
                continue

            df2 = df.append(healthy_df, ignore_index=True)
            fpr,tpr,_ = roc_curve(df2['status'],df2[i])
            AUC = auc(fpr,tpr)
            current_dict[group] = AUC
            del(AUC,group,fpr,tpr)
            
        for group,df in cancer_df.groupby('Stage'):
            if group == '0' or group == 'X':
                continue
            df2 = df.append(healthy_df, ignore_index=True)
            fpr,tpr,_ = roc_curve(df2['status'],df2[i])
            AUC = auc(fpr,tpr)
            current_dict[group] = AUC
            del(AUC,group,fpr,tpr)
            
        for group,df in cancer_df.groupby('tfx_group'):
            df2 = df.append(healthy_df, ignore_index=True)
            fpr,tpr,_ = roc_curve(df2['status'],df2[i])
            AUC = auc(fpr,tpr)
            current_dict[group] = AUC
            del(AUC,group,fpr,tpr)
            
        AUCs = AUCs.append(pd.Series(current_dict), ignore_index=True)
        
    CIs = pd.DataFrame([AUCs.median(), AUCs.quantile(.025), AUCs.quantile(.975)]).T
    CIs = CIs.rename(columns = {'Unnamed 0':'median'})    
    return(AUCs,CIs)


# In[ ]:


print('analysis name',args.name)
print('importing data')
data,features  = import_data(args.in_file)
sys.stdout.flush()


# In[ ]:


print('sample type:')
print(data['sample_type'].value_counts(dropna = False))
print('\n')

print('stages:')
print(data[data['status']==1]['Stage'].value_counts(dropna = False))
print('\n')

print('tfx:')
print(data['tfx_group'].value_counts(dropna = False))


# In[ ]:


print('running '+str(args.iterations)+' logreg bootstrap iterations')
hyperparameters = {'C': [0.0001, 0.001,0.01,0.1,1,10,100,1000]}

probabilities,c_vals,coefs,num_pcs,train_indexes = run_bootstrap_with_PCA(data,args.iterations,features,args.report_interval,hyperparameters)    
sys.stdout.flush()


# In[ ]:


print('Getting AUC')
sys.stdout.flush()
AUCs,CIs = get_AUC(probabilities,data,args.iterations)


# In[ ]:


probabilities.to_csv(args.out_dir+'/'+args.name+'.probabilities.txt', sep='\t', float_format='%.5f')
pd.Series(c_vals).to_csv(args.out_dir+'/'+args.name+'.c_values.txt', sep='\t', header = False, index=False)
coefs.to_csv(args.out_dir+'/'+args.name+'.coefs.txt', sep='\t', float_format='%.5f')
pd.Series(num_pcs).to_csv(args.out_dir+'/'+args.name+'.num_pcs.'+str(args.fraction_variance)+'.txt', sep='\t', header = False, index=False)
pd.Series(train_indexes).to_csv(args.out_dir+'/'+args.name+'.train_indexes.txt', sep='\t', header = False, index=False)

AUCs.to_csv(args.out_dir+'/'+args.name+'.AUC.txt', sep='\t', index = False, float_format='%.5f')
CIs.to_csv(args.out_dir+'/'+args.name+'.CI.txt', sep='\t', float_format = '%.5f')


# In[ ]:


bins = list(np.log10(hyperparameters['C']))+[np.log10(max(hyperparameters['C'])*10)]
plt.hist(np.log10(c_vals),bins = bins, label=hyperparameters['C'],  align = 'left', rwidth = .8);
plt.xlabel('log10 C value')
plt.ylabel('Frequency')
plt.savefig(args.out_dir+'/'+args.name+'.cvals.pdf')   
print('done')


# In[ ]:


plt.hist(num_pcs)
plt.ylabel('frequency')
plt.xlabel('number of pcs')
plt.savefig(args.out_dir+'/'+args.name+'.pcs.'+str(args.fraction_variance)+'.pdf') 


# In[ ]:





# In[ ]:





# In[ ]:




