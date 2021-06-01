#view the results
#Usage: python results.py resultfile

import pandas as pd
import sys
from sklearn.metrics import accuracy_score, roc_auc_score as roc, auc, classification_report,balanced_accuracy_score, matthews_corrcoef as mcc
df = pd.read_csv(sys.argv[1], header=None, sep='\t',)
y_true = df[1]
y_pred = df[2]
y_predprob = df[3]
#print(y_test)
print("Accuracy:",accuracy_score(y_true, y_pred))
print("Balanced accuracy:",balanced_accuracy_score(y_true, y_pred))
#print(auc(y_true, y_predprob))
print(classification_report(y_true, y_pred))
print("MCC:",mcc(y_true, y_pred))
print("ROC_AUC:",roc(y_true, y_predprob))
