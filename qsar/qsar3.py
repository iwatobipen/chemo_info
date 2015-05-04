#! /usr/env/python

import numpy as np
import pandas as pd
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import PandasTools
from rdkit.ML.Descriptors import MoleculeDescriptors

from sklearn.cross_decomposition import PLSCanonical, PLSRegression

from sklearn import cross_validation
from sklearn import metrics

trainset = PandasTools.LoadSDF("solubility.train.sdf")
testset =  PandasTools.LoadSDF("solubility.test.sdf")

# calclator
nms = [ x[0] for x in Descriptors._descList ]

def calc_desc( mol ):
    calc = MoleculeDescriptors.MolecularDescriptorCalculator( nms )
    res = calc.CalcDescriptors( mol )
    return res

trainset["desc"]=trainset.ROMol.apply(calc_desc)
train_desc = pd.DataFrame([ i for i in trainset.desc ])
print train_desc.shape

