import sys, cPickle
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

from sklearn.ensemble import RandomForestClassifier
from sklearn import svm
from sklearn.naive_bayes import GaussianNB
from sklearn import neighbors
from sklearn import cross_validation
from sklearn import tree, metrics
from sklearn.cross_validation import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from sklearn.grid_search import GridSearchCV
from sklearn.metrics import classification_report

trainset_file = sys.argv[1]
testset_file = sys.argv[2]
trainset = [ x for x in Chem.SDMolSupplier(trainset_file) ]
testset = [ x for x in Chem.SDMolSupplier(testset_file) ]

classes = { '(A) low':0, '(B) medium':1, '(C) high':2 }
desc_names = [ x[0] for x in Descriptors._descList ]

def calclator( mols ):
    nms = [ x[0] for x in Descriptors._descList ]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator( nms )
    desc = [ calc.CalcDescriptors(mol) for mol in mols ]
    return desc

#################################
# make training set and test set#
#################################
trainActs = np.array([ classes[x.GetProp( 'SOL_classification' )] for x in trainset ])
testActs = np.array([ classes[x.GetProp( 'SOL_classification' )] for x in testset ])
trainDescrs = calclator( trainset )
testDescrs = calclator( testset )
cPickle.dump(trainDescrs,file('trainDescrs.pkl','wb+'))
cPickle.dump(testDescrs,file('testDescrs.pkl','wb+'))


print "RandomForest"
clf_RF = RandomForestClassifier( n_estimators=100, max_depth=5,  min_samples_split=5, random_state=0,n_jobs=1 )
scores_RF = cross_validation.cross_val_score(clf_RF, trainDescrs, trainActs)
print scores_RF
nclf_RF = RandomForestClassifier( n_estimators=100, max_depth=5,  min_samples_split=5, random_state=0,n_jobs=1 )
nclf_RF = nclf_RF.fit( testDescrs, testActs )
pred_RF = nclf_RF.predict( testDescrs )
print metrics.confusion_matrix( testActs, pred_RF )
accuracy = nclf_RF.score( testDescrs, testActs )
print accuracy
print ""
print ""

print "Naive_Bayse"
gnb = GaussianNB()
clf_NB = gnb.fit( trainDescrs, trainActs )
pred_NB = clf_NB.predict( testDescrs )
print metrics.confusion_matrix( testActs, pred_NB )
accuracy = clf_NB.score( testDescrs, testActs )
print accuracy

print ""
print ""

print "K-NN"
clf_KNN = neighbors.KNeighborsClassifier( 15, weights="uniform" )
clf_KNN.fit( trainDescrs, trainActs )
pred_KNN = clf_KNN.predict( testDescrs )
print metrics.confusion_matrix( testActs, pred_KNN )
accuracy = clf_KNN.score( testDescrs, testActs )
print accuracy

print ""
print ""


print "DecsionTreeClassifier"
clf_DT = tree.DecisionTreeClassifier()
clf_DT.fit( trainDescrs, trainActs )
pred_DT = clf_DT.predict( testDescrs )
print metrics.confusion_matrix( testActs, pred_DT )
accuracy = clf_DT.score( testDescrs, testActs )
print accuracy
export_file = tree.export_graphviz( clf_DT, out_file='tree.dot', feature_names = desc_names)

print ""
print ""

print "SVM"
clf_svm = svm.SVC( gamma=0.001, C=100. )
clf_svm = clf_svm.fit( trainDescrs, trainActs )
pred_SVM = clf_svm.predict( testDescrs )
print metrics.confusion_matrix( testActs, pred_SVM )
accuracy = clf_svm.score( testDescrs, testActs )
print accuracy

print "finished"
