import pandas as pd
import sys
from sklearn.metrics import accuracy_score, roc_auc_score as roc, auc, classification_report,balanced_accuracy_score, matthews_corrcoef as mcc
df = pd.read_csv(sys.argv[1], header=0)
for col in df.columns:
    print(col)
y_true = df['target']
y_pred = df['pred']
y_predprob = df['pred-score']
#print(y_test)
print(accuracy_score(y_true, y_pred))
print(balanced_accuracy_score(y_true, y_pred))
print(auc(y_true, y_predprob))
print(classification_report(y_true, y_pred))
print(mcc(y_true, y_pred))
print(roc(y_true, y_predprob))
