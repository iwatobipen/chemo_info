
import nbviz
import numpy as np
import sys
import maccskey
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from sklearn.naive_bayes import BernoulliNB


def calc_MACCS_fp( mol ):
	mol_fp =list( MACCSkeys.GenMACCSKeys( mol ).GetOnBits() )
	mol_fp_vec = np.zeros( 167, )
	mol_fp_vec[ mol_fp ] = 1
	return mol_fp_vec

def make_fp_array( mols ):
	fp_array = [ calc_MACCS_fp( mol ) for mol in mols ]
	return fp_array





mols = [ mol for mol in Chem.SDMolSupplier( sys.argv[1] ) ]
X = make_fp_array( mols )
Y = [ float(mol.GetProp( "Class" )) for mol in mols ]

model = BernoulliNB( alpha=0.1 )
model.fit( X[1:], Y[1:] )
conditional_probs = np.exp( model.feature_log_prob_ )
prior = np.exp( model.class_log_prior_[1] )
print 'condtional feature prob', conditional_probs
print 'class prior', prior
nbviz.visualize_model( conditional_probs, prior, names=maccskey.names, groups=maccskey.groups )

nbviz.visualize_prediction( X[0], conditional_probs, prior, names=maccskey.names, groups=maccskey.groups )

