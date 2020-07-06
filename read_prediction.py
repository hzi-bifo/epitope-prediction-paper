import sys
sys.path.append('./')
from utility.file_utility import FileUtility
import numpy as np
path2res = sys.argv[1]
files = sys.argv[1]
import warnings
warnings.filterwarnings("ignore")


# In[8]:


def get_cv_res(filename):
    [label_set, conf, best_score_, best_estimator_, cv_results_,
        best_params_, pred] = FileUtility.load_obj(filename)
    res = dict()
    print (conf)
    #print (cv_results_.keys())
    idx = np.argmax(cv_results_['mean_test_f1_macro'])
    res['f1_macro'] = np.round(cv_results_['mean_test_f1_macro'][idx], 2)
    res['f1_macro*'] = str(np.round(cv_results_['mean_test_f1_macro'][idx], 2)) + \
        ' $\pm$ ' + str(np.round(cv_results_['std_test_f1_macro'][idx], 2))
    res['f1_micro'] = str(np.round(cv_results_['mean_test_f1_micro'][idx], 2)) + \
        ' $\pm$ ' + str(np.round(cv_results_['std_test_f1_micro'][idx], 2))
    res['precision_micro'] = str(np.round(cv_results_['mean_test_precision_micro'][idx], 2)) + \
        ' $\pm$ ' + \
        str(np.round(cv_results_['std_test_precision_micro'][idx], 2))
    res['precision_macro'] = str(np.round(cv_results_['mean_test_precision_macro'][idx], 2)) + \
        ' $\pm$ ' + \
        str(np.round(cv_results_['std_test_precision_macro'][idx], 2))
    res['recall_micro'] = str(np.round(cv_results_['mean_test_recall_micro'][idx], 2)) + \
        ' $\pm$ ' + str(np.round(cv_results_['std_test_recall_micro'][idx], 2))
    res['recall_macro'] = str(np.round(cv_results_['mean_test_recall_macro'][idx], 2)) + \
        ' $\pm$ ' + str(np.round(cv_results_['std_test_recall_macro'][idx], 2))
    #res['accuracy']=str(np.round(cv_results_['mean_test_accuracy'][idx],2))+ ' $\pm$ ' + str(np.round(cv_results_['std_test_accuracy'][idx],2))
    res['file'] = file
    res['auc_macro'] = str(conf['auc_macro'])
    res['score'] = str(best_score_)
    return res


res = dict()
for file in files:
    k = ''.join(file.split('/')[-1].split('_')[0:2])
    s = (file.split('/')[-1].split('_')[2])
    if k not in res:
        res[k] = dict()
        res[k][s] = get_cv_res(file)


# In[9]:


keys = list(res.keys())
print (keys)
keys.sort()
for k in keys:

    try:
        keys2 = [int(x) for x in list(res[k].keys())]
    except:
        keys2 = list(res[k].keys())
    keys2.sort()
    max_val = -1
    max_arg = -1
    for k2 in keys2:
        k2 = str(k2)
        if res[k][k2]['f1_macro'] > max_val:
            max_arg = k2
            max_val = res[k][k2]['f1_macro']
    print (' & '.join([str(k), str(max_arg), res[k][k2]['auc_macro'], res[k][k2]['score'], res[k][k2]['precision_micro'], res[k][k2]['recall_micro'],
                       res[k][k2]['f1_micro'], res[k][k2]['precision_macro'], res[k][k2]['recall_macro'], res[k][k2]['f1_macro*']]) + '\\ \\hline')
