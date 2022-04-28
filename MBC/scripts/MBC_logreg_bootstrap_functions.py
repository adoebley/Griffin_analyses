#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def import_data(in_file, status_col, min_tfx, min_cov):
    import pandas as pd
    import numpy as np
    from sklearn.preprocessing import StandardScaler
    
    data = pd.read_csv(in_file, sep='\t')
    data = data.set_index('sample')

    #get features and exclude all other columns
    features = data.columns[(data.columns.str.startswith('central_cov')) | (data.columns.str.startswith('mean_cov')) | (data.columns.str.startswith('amplitude'))]
    print('Features',len(features))

    #get status
    data['status'] = np.where(data[status_col]=='+',1,0)

    #filter data
    data = data[(data['tumor_fraction']>=min_tfx) & (data['ulp_wgs_coverage']>=min_cov)]
    data = data.sort_index()

    print('Total samples:',len(data))
    
#     display(pd.DataFrame(data[status_col].value_counts()))

#     high_tfx = data[(data['tumor_fraction']>=0.1) & (data['ulp_wgs_coverage']>=min_cov)]
#     print('High tfx samples:',len(high_tfx))
#     display(pd.DataFrame(high_tfx[status_col].value_counts()))
#     del(high_tfx)

#     low_tfx = data[(data['tumor_fraction']<0.1) & (data['ulp_wgs_coverage']>=min_cov)]
#     print('Low tfx samples:',len(low_tfx))
#     display(pd.DataFrame(low_tfx[status_col].value_counts()))
#     del(low_tfx)

    #scale data
    scaler = StandardScaler()
    scaler.fit(data[features])
    data[features] = scaler.transform(data[features])
    data[features].mean()
    
    return(data,features)


# In[ ]:


#bootstrapping
def run_bootstrap(data,iterations,features,report_interval):
    import time
    import sys
    import pandas as pd
    import numpy as np
    from sklearn.model_selection import StratifiedKFold
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GridSearchCV
    from matplotlib import pyplot as plt
    
    hyperparameters = {'C': [0.00001, 0.0001, 0.001,0.01,0.1,1,10,100]}

    start_time = time.time()

    probabilities = pd.DataFrame(index=data.index)
    c_vals = []
    coefs = pd.DataFrame(index=features)

    #loop for each iteration
    countup = 0
    for i in range(iterations):
        if i%report_interval==0:
            print(i, time.time()-start_time)
            #prevent dfs from becoming too fragmented
            probabilities = probabilities.copy()
            coefs = coefs.copy()   
            sys.stdout.flush() 

        patients = pd.Series(data['patient_id'].unique())

        good_split = 0
        while good_split == 0:
            #bootstrap a training set with replacement
            training_ids = patients.sample(len(patients), replace = True, random_state = 100+countup)

            #get bootstrapped training set, if a patient ID is included in the training_ids set j times, include all samples from that patient j times
            training = pd.DataFrame()
            #group the training patient IDs by number of times they are observed in the bootstrapped training_ids
            for j,df in pd.DataFrame(training_ids.value_counts().rename('count')).groupby(by = 'count'):

                #identify the samples from this group of patients
                current_data = data[data['patient_id'].isin(df.index)]

                #copy the training samples so that they appear j times in training dataframe
                current_training = pd.DataFrame()
                for k in range(j):
                    current_training = current_training.append(current_data)
                training = training.append(current_training)


            #the test set is all samples that aren't seen in the training data
            test = data[~(data.index.isin(training.index))]

            #print(len(training),len(training.index.unique()))
            #print(len(test),len(test.index.unique()))

            #check to make sure first time point low tumor fraction samples for both classes are included in the test set
            if len(test[(test['tumor_fraction']<.1) & (test['first_passing_sample']==1)]['status'].unique())!=2:
                print('Skipping',i)
                countup += 1
            else:
                good_split = 1

        #countup will get ahead of i if it has to skip bad train-test splits
        countup +=1

        #10 fold cross validation on the training set
        cv = StratifiedKFold(n_splits=10, shuffle=True, random_state =100+countup) 

        model = LogisticRegression(class_weight='balanced', max_iter=500)
        search = GridSearchCV(estimator=model, param_grid=hyperparameters, cv=cv, n_jobs = 1)
        search.fit(training[features], training['status'])
        best_C = search.best_params_['C']

        #train a new model on the full training dataset (is this the same as refit...?)
        model = LogisticRegression(class_weight='balanced', max_iter=500, C=best_C)
        model.fit(training[features], training['status'])

        #predict the test data
        pred = model.predict(test[features])
        prob = model.predict_proba(test[features])


        #collect metrics
        current_output = pd.DataFrame(test[['status']])#.reset_index()
        current_output['probability']=prob[:,1]
        current_output['prediction']=pred
        current_output['accuracy'] = np.where(current_output['prediction']==current_output['status'],1,0)
  
        #save results
        probabilities[i] = current_output['probability']
        c_vals.append(best_C)
        coefs[i] = pd.Series(model.coef_[0], index = features)

    probabilities = probabilities.merge(data[['status']], left_index=True, right_index=True)
    
    return(probabilities,c_vals,coefs)


# In[ ]:


# def get_accuracy(probabilities,data,iterations):
#     import numpy as np
#     probabilities = probabilities.merge(data[['first_passing_sample','tumor_fraction']], left_index=True, right_index=True)
    
#     #high tfx per patient accuracy
#     accuracy = []
#     for i in range(iterations):
#         current = probabilities[~(probabilities[i].isnull())][[i,'status','tumor_fraction','first_passing_sample']]
#         current = current[(current['tumor_fraction']>=0.1) & (current['first_passing_sample']==1)]
#         accuracy.append(sum(np.round(current[i])==current['status'])/len(current))

#     high_tfx_accuracy = np.mean(accuracy)
#     print('high tfx per patient:',np.round(np.mean(accuracy),3))
    
#     #low tfx per patient accuracy
#     accuracy = []
#     for i in range(iterations):
#         current = probabilities[~(probabilities[i].isnull())][[i,'status','tumor_fraction','first_passing_sample']]
#         current = current[(current['tumor_fraction']<0.1) & (current['first_passing_sample']==1)]
#         accuracy.append(sum(np.round(current[i])==current['status'])/len(current))
    
#     low_tfx_accuracy = np.mean(accuracy)
#     print('low tfx per patient:',np.round(np.mean(accuracy),3))
    
#     #per patient accuracy for all samples
#     accuracy = []
#     for i in range(iterations):
#         current = probabilities[~(probabilities[i].isnull())][[i,'first_passing_sample','status']]
#         current = current[(current['first_passing_sample']==1)]
#         accuracy.append(sum(np.round(current[i])==current['status'])/len(current))

#     all_tfx_accuracy = np.mean(accuracy)
#     print('all tfx per patient:',np.round(all_tfx_accuracy,3))
    
#     return(high_tfx_accuracy,low_tfx_accuracy,all_tfx_accuracy)


# In[3]:


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


# In[ ]:




