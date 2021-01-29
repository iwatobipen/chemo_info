#
import numpy as np
from matplotlib import cm
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Draw import SimilarityMaps
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
import seaborn as sns
import os, sys
dataset = open( "data_cdk.txt", "r" )
# remove header
dataset = [ line.rstrip().split("\t") for line in dataset ][1:]
mols = [ Chem.MolFromSmiles( line[1] ) for line in dataset ]
Y = []
# label active / non active
for line in dataset:
    if float( line[2] ) > 7.0:
        Y.append( 1 )
    else:
        Y.append( -1 )
Y = np.asarray( Y )

trainmols, testmols, trainy, testy = train_test_split( mols, Y, test_size = 0.1, random_state=90  )
# calc fingerprint
trainfps = [ AllChem.GetMorganFingerprintAsBitVect( mol,2 ) for mol in trainmols ]
testfps =  [ AllChem.GetMorganFingerprintAsBitVect( mol,2 ) for mol in testmols ]

trainx = []
for fp in trainfps:
    arr = np.zeros( (1,) )
    DataStructs.ConvertToNumpyArray( fp, arr )
    trainx.append( arr )

testx = []
for fp in testfps:
    arr = np.zeros( (1,) )
    DataStructs.ConvertToNumpyArray( fp, arr )
    testx.append( arr )

#print(len(mols),len(trainx), len(trainy))
cls = SVC( probability=True, C=100 )
cls.fit( trainx, trainy )
def getProba( fp, probabilityFunc ):
    return probabilityFunc([fp])[0][1]

def mapperfunc( mol ):
    fig, weight = SimilarityMaps.GetSimilarityMapForModel( mol, SimilarityMaps.GetMorganFingerprint, lambda x: getProba( x, cls.predict_proba), colorMap=cm.bwr  )
    fp = AllChem.GetMorganFingerprintAsBitVect( mol, 2 )
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray( fp, arr )
    res = cls.predict( [arr] )
    smi = Chem.MolToSmiles( mol )
    print(res[0])
    if res[0] == 1:
        fig.savefig( "res/act_"+smi+"_.png", bbox_inches = "tight" )
    else:
        fig.savefig("res/nonact_"+smi+"_.png", bbox_inches = "tight" )

for mol in testmols:
    mapperfunc( mol )

res= cls.predict(testx)
from sklearn.metrics import confusion_matrix
print(confusion_matrix(testy, res))
