from Bio import SeqIO
from pydpi.pypro import PyPro
from make_representations.sequencelist_representation import SequenceKmerRep, SequenceKmerEmbRep
from sklearn import svm, preprocessing
import sys
import numpy as np
import warnings
warnings.filterwarnings("ignore")
from sklearn import svm, datasets, metrics
from sklearn import preprocessing
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score, roc_curve, precision_recall_curve, precision_recall_fscore_support, roc_auc_score
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import classification_report
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold, GridSearchCV, StratifiedKFold
# from sklearn.cross_validation import 
import numpy as np
import sys
import pylab as pl
import matplotlib.pyplot as plt
#from sklearn.grid_search import GridSearchCV
stdsc = StandardScaler()
from matplotlib.colors import ListedColormap
from sklearn.neural_network import MLPClassifier
from scipy import interp
from sklearn.metrics import matthews_corrcoef as mcc
from pydpi.pypro import PyPro
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
# from classifier.classical_classifiers import RFClassifier, SVM
# from make_representations.sequencelist_representation import SequenceKmerRep, SequenceKmerEmbRep
from sklearn.feature_selection import VarianceThreshold
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import precision_score, recall_score
import argparse
import pickle

protein = PyPro()


def readAAP(file):  #read AAP features from the AAP textfile
    try:
        aapdic = {}
        aapdata = open(file, 'r')
        for l in aapdata.readlines():
            aapdic[l.split()[0]] = float(l.split()[1])
        aapdata.close()
        return aapdic
    except:
        print("Error in reading AAP feature file. Please make sure that the AAP file is correctly formatted")
        sys.exit()


def readAAT(file):  #read AAT features from the AAT textfile
    try:
        aatdic = {}
        aatdata = open(file, 'r')
        for l in aatdata.readlines():
            aatdic[l.split()[0][0:3]] = float(l.split()[1])
        aatdata.close()
        return aatdic
    except:
        print("Error in reading AAT feature file. Please make sure that the AAT file is correctly formatted")
        sys.exit()


def aap(pep, aapdic, avg):  #return AAP features for the peptides
    feature=[]
    for a in pep:
        #print(a)
        if int(avg) == 0:
            score = []
            count = 0
            for i in range(0, len(a) - 1):
                try:
                    score.append(round(float(aapdic[a[i:i + 2]]), 4))
                    # score += float(aapdic[a[i:i + 3]])
                    count += 1
                except KeyError:
                    # print(a[i:i + 3])
                    score.append(float(-1))
                    # score += -1
                    count += 1
                    continue
            # averagescore = score / count
            feature.append(score)
        if int(avg) == 1:
            score = 0
            count = 0
            for i in range(0, len(a) - 1):
                try:
                    score += float(aapdic[a[i:i + 2]])
                    count += 1
                except KeyError:
                    score += -1
                    count += 1
                    continue
            if count != 0:
                averagescore = score / count
            else:
                averagescore = 0
            feature.append(round(float(averagescore), 4))
    return feature


def aat(pep, aatdic, avg):  #return AAT features for the peptides
    feature = []
    for a in pep:
        if int(avg) == 0:
            # print(a)
            score = []
            count = 0
            for i in range(0, len(a) - 2):
                try:
                    score.append(round(float(aatdic[a[i:i + 3]]), 4))
                    # score += float(aapdic[a[i:i + 3]])
                    count += 1
                except KeyError:
                    # print(a[i:i + 3])
                    score.append(float(-1))
                    # score += -1
                    count += 1
                    continue
            # averagescore = score / count
            feature.append(score)
        if int(avg) == 1:
            score = 0
            count = 0
            for i in range(0, len(a) - 2):
                try:
                    score += float(aatdic[a[i:i + 3]])
                    count += 1
                except KeyError:
                    score += -1
                    count += 1
                    continue
            # print(a, score)
            if count != 0:
                averagescore = score / count
            else:
                averagescore = 0
            feature.append(round(float(averagescore), 4))
    return feature


def CTD(pep):  #Chain-Transition-Ditribution feature
    feature = []
    for seq in pep:
        protein.ReadProteinSequence(seq)
        ctd = protein.GetCTD()
        feature.append(list(ctd.values()))
    return feature


def AAC(pep): # Single Amino Acid Composition feature
    feature = []
    for seq in pep:
        protein.ReadProteinSequence(seq)
        aac = protein.GetAAComp()
        feature.append(aac)
    return feature


def DPC(pep): # Dipeptide Composition feature
    feature = []
    for seq in pep:
        protein.ReadProteinSequence(seq)
        dpc = protein.GetDPComp()
        feature.append(list(dpc.values()))
    return feature


def kmer(pep, k): # Calculate k-mer feature
    feature = SequenceKmerRep(pep, 'protein', k)
    return feature


def protvec(pep, k, file): #Calculate ProtVec representation
    feature = SequenceKmerEmbRep(file, pep, 'protein', k)
    return feature


def readseq(file):  #read the sequence from the fasta file
    try:
        sequence = SeqIO.read(file, "fasta")
        for i in sequence.seq:
            #print(i)
            if i in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'] :
                continue
            else:
                print("Invalid amino acid code found. Please enter sequences with only 20 aa code.")
                sys.exit()
        return(str(sequence.seq))
    except ValueError:
        print("Please enter a valid fasta file")
        sys.exit()


def peptides(seq):  #return peptides of length 20 from the sequence
    pep = []
    i=0
    while i <= len(seq):
        if i+20 > len(seq):
            pep.append(seq[i:len(seq)])
        else:
            pep.append(seq[i:i+20])
        i = i + 20
    #print(pep)
    return pep


def readpeptides(posfile, negfile): # return the peptides from input peptide list file 
    posdata = open(posfile,'r')
    pos = []
    for l in posdata.readlines():
        if l[0] == '#':
            continue
        else:
            pos.append(l.strip)
    negdata = open(negfile,'r')
    neg = []
    for l in negdata.readlines():
        if l[0] == '#':
            continue
        else:
            neg.append(l.strip)
    return pos, neg


def combinefeature(pep):
    aapdic = readAAP("./aap-general.txt.normal")
    aatdic = readAAT("./aat-general.txt.normal")
    f_aap = np.array(aap(pep, aapdic, 1))
    #print(f_aap)
    f_aat = np.array(aat(pep, aatdic, 1))
    #print(f_aat)
    f_aac = np.array(AAC(pep))
    #print(f_aac)
    f_kmer = np.array(kmer(pep, 4).X.toarray())
    #print(f_kmer)
    f_protvec = np.array(protvec(pep, 4, '../epitope-prediction-master/protvec/sp_sequences_4mers_vec.txt').embeddingX)
    #print(f_protvec)
    return np.column_stack((f_aat,f_aac,f_kmer,f_protvec))


def gridsearch(x, y, cv):
    grid_search = GridSearchCV(SVC(kernel='rbf', probability=True),
                               param_grid={'C': [1000, 500, 200, 100, 50,
                                                 20, 10, 2, 1, 0.2, 0.5,
                                                 0.01, 0.02, 0.05, 0.001],
                                           'gamma': [1000, 500, 200, 100,
                                                     50, 20, 10, 5, 2, 1,
                                                     0.2, 0.5, 0.01, 0.02,
                                                     0.05, 0.001, 0.0001]},
                               scoring='accuracy', cv=cv, n_jobs=40)
    grid_search.fit(x, y)
    return grid_search


def train(peptides, features, target):
    scaling = StandardScaler()
    scaling.fit(features)
    x = scaling.transform(features)
    y = target
    cv = StratifiedKFold(n_splits = 5)
    model = gridsearch(x, y, cv)
    pickle.dump(scaling, open("scaling.pickle", "wb"))
    pickle.dump(model, open("model.pickle", "wb") 
    print("Best parameters: ", grid_search.best_params_)
    print("Best accuracy: :", grid_search.best_score_)
    print("Best 

    


def readmodel(mlfile):
    try:
        return pickle.load(open(mlfile, 'rb'))
    except:
        print("Error in reading model file")
        sys.exit()


def predict(model, features):
    try:
        return model.predict_proba(features)
    except:
        print("Error in predicting epitopes.")
        sys.exit()





def scoremodel(file, mlfile):
    sequence = readseq(file)
    pep = peptides(sequence)
    features = combinefeature(pep)
    #print(len(features[0]))
    model = readmodel(mlfile)
    return pep, predict(model, features)


if __name__ == "__main__":
    peptide_list, pred_probability = scoremodel("./input/example.fasta", "./model/svm-model.pickle" )
    print("List of predicted epitopes:")
    for i in range(len(pred_probability)):
        if pred_probability[i][1] >= 0.5:
            print(peptide_list[i], pred_probability[i][1])
