#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import sys
import os
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np


# In[2]:


print('starting')
sys.stdout.flush()


# In[3]:


# %matplotlib inline

# cmd = ['--in_file','../number_of_sites_analysis/merged_data/30000_reformatted.txt',
#        '--name','30000',
#         '--status_column','revisions_ER_status_binary',
#       '--iterations','20',
#       '--out_dir','tmp/tmp/',
#       '--report_interval','50',
#        '--fraction_variance','.8']


# In[4]:


parser = argparse.ArgumentParser()

parser.add_argument('--in_file', help='/path/to/data_reformatted.txt', required=True)
parser.add_argument('--name', help='name for outputs', required=True)
parser.add_argument('--status_column', help='column containing the status', required=True)
parser.add_argument('--iterations', help='number of iterations', required=True, type = int)
parser.add_argument('--out_dir', help='output directory', required=True)
parser.add_argument('--report_interval', help='how frequently to print progress', default = 50, type = int)
parser.add_argument('--fraction_variance', help='what fraction of the variance should be explained by the PCs used in the model', type = float)

args = parser.parse_args()


# In[5]:


#create directory to export results
if not os.path.exists(args.out_dir):
    os.mkdir(args.out_dir)


# In[6]:


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
    data['status'] = data[args.status_column]
    data['status'] = data['status'].replace('+',1).replace('-',0)
    
    return(data,features)


# In[7]:


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


# In[8]:


def get_accuracy_AUC(probabilities,data,iterations):
    #get AUC and accuracy for each bootstrap
    from sklearn.metrics import roc_curve,auc
    import pandas as pd
    import numpy as np

    AUCs = pd.DataFrame()
    accuracies = pd.DataFrame()

    probabilities = probabilities.merge(data[['first_passing_sample','tumor_fraction']], left_index=True, right_index=True)
    probabilities = probabilities[probabilities['first_passing_sample']==1]

    for i in range(iterations):
        current_accuracy_dict = {}
        current_AUC_dict = {}
        
        current = probabilities[~(probabilities[i].isnull())][['status','tumor_fraction',i]].copy()

        current['accuracy'] = np.where(np.round(current[i])==current['status'],1,0)
        
        low_tfx = current[(current['tumor_fraction']<0.1)].copy()
        high_tfx = current[(current['tumor_fraction']>=0.1)].copy()

        for group,df in zip(['overall','high_tfx','low_tfx'],[current,high_tfx,low_tfx]):
            fpr,tpr,_ = roc_curve(df['status'],df[i])

            AUC = auc(fpr,tpr)
            current_AUC_dict[group] = AUC

            accuracy = sum(np.round(df[i])==df['status'])/len(df)
            current_accuracy_dict[group] = accuracy
            
        accuracies = accuracies.append(pd.Series(current_accuracy_dict), ignore_index=True)
        AUCs = AUCs.append(pd.Series(current_AUC_dict), ignore_index=True)
        

    AUC_CIs = pd.DataFrame([AUCs.median(), AUCs.quantile(.025), AUCs.quantile(.975)]).T
    AUC_CIs = AUC_CIs.rename(columns = {'Unnamed 0':'median'})    

    accuracy_CIs = pd.DataFrame([accuracies.median(), accuracies.quantile(.025), accuracies.quantile(.975)]).T
    accuracy_CIs = accuracy_CIs.rename(columns = {'Unnamed 0':'median'})    

    return(accuracies,accuracy_CIs,AUCs,AUC_CIs)


# In[9]:


print('analysis name',args.name)
print('importing data')
data,features  = import_data(args.in_file)
sys.stdout.flush()


# In[10]:


print('status')
print(data['status'].value_counts())
    
print('tfx:')
low_tfx = data[(data['tumor_fraction']>=0.05) & (data['tumor_fraction']<0.1)]
high_tfx = data[(data['tumor_fraction']>=0.1)]

print('low tfx:',len(low_tfx))
print('high tfx:',len(high_tfx))

del(high_tfx,low_tfx)


# In[11]:


print('running '+str(args.iterations)+' logreg bootstrap iterations')
hyperparameters = {'C': [0.0001, 0.001,0.01,0.1,1,10,100,1000]}

probabilities,c_vals,coefs,num_pcs,train_indexes = run_bootstrap_with_PCA(data,args.iterations,features,args.report_interval,hyperparameters)    
sys.stdout.flush()


# In[12]:


print('Getting AUC')
sys.stdout.flush()
accuracies,accuracy_CIs,AUCs,AUC_CIs = get_accuracy_AUC(probabilities,data,args.iterations)


# In[13]:


probabilities.to_csv(args.out_dir+'/'+args.name+'.probabilities.txt', sep='\t', float_format='%.5f')
pd.Series(c_vals).to_csv(args.out_dir+'/'+args.name+'.c_values.txt', sep='\t', header = False, index=False)
coefs.to_csv(args.out_dir+'/'+args.name+'.coefs.txt', sep='\t', float_format='%.5f')
pd.Series(num_pcs).to_csv(args.out_dir+'/'+args.name+'.num_pcs.'+str(args.fraction_variance)+'.txt', sep='\t', header = False, index=False)
pd.Series(train_indexes).to_csv(args.out_dir+'/'+args.name+'.train_indexes.txt', sep='\t', header = False, index=False)

accuracies.to_csv(args.out_dir+'/'+args.name+'.accuracy.txt', sep='\t', index = False, float_format='%.5f')
accuracy_CIs.to_csv(args.out_dir+'/'+args.name+'.accuracy_CI.txt', sep='\t', float_format = '%.5f')
AUCs.to_csv(args.out_dir+'/'+args.name+'.AUC.txt', sep='\t', index = False, float_format='%.5f')
AUC_CIs.to_csv(args.out_dir+'/'+args.name+'.AUC_CI.txt', sep='\t', float_format = '%.5f')


# In[14]:


bins = list(np.log10(hyperparameters['C']))+[np.log10(max(hyperparameters['C'])*10)]
plt.hist(np.log10(c_vals),bins = bins, label=hyperparameters['C'],  align = 'left', rwidth = .8);
plt.xlabel('log10 C value')
plt.ylabel('Frequency')
plt.savefig(args.out_dir+'/'+args.name+'.cvals.pdf')   
print('done')


# In[32]:


plt.hist(num_pcs, range = (min(num_pcs),max(num_pcs)+1), bins = max(num_pcs) - min(num_pcs)+1, align = 'left', rwidth=.8)
plt.ylabel('frequency')
plt.xlabel('number of pcs')
plt.savefig(args.out_dir+'/'+args.name+'.pcs.'+str(args.fraction_variance)+'.pdf') 


# In[33]:


print('AUCs')
print(AUC_CIs)


# In[35]:


print('Accuracies')
print(accuracy_CIs)


# In[ ]:





# In[ ]:




