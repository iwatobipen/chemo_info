#! /usr/bin/python

from sklearn.cross_decomposition import PLSCanonical, PLSRegression
from sklearn import metrics
import numpy as np

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
nms = [ x[0] for x in Descriptors._descList ]
def calculator( mols ):
    calc = MoleculeDescriptors.MolecularDescriptorCalculator( nms )
    res = [ calc.CalcDescriptors( mol ) for mol in mols ]
    return res

trainMols = [ mol for mol in Chem.SDMolSupplier("solubility.train.sdf") ]
testMols =  [ mol for mol in Chem.SDMolSupplier("solubility.test.sdf") ]

trainDescrs = calculator( trainMols )
testDescrs = calculator( testMols )

trainActs = np.array([ float( mol.GetProp('SOL') ) for mol in trainMols  ])
testActs = np.array([ float( mol.GetProp('SOL') ) for mol in testMols  ])

pls2 = PLSRegression( n_components = 15 )
pls2.fit( trainDescrs, trainActs )

sol_pred = pls2.predict( testDescrs )
print type(sol_pred)
print type(trainActs)
print metrics.r2_score( testActs, sol_pred )

"""
for i in range(len(sol_pred)):
    print testActs[i], sol_pred[i]

"""
