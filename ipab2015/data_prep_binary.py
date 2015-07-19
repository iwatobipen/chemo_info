from sklearn.cross_validation import train_test_split

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs

import gzip
import cPickle
import numpy as np

mols = [ mol for mol in Chem.SDMolSupplier( 'dataset/dataset_sdf_clean.sdf' ) ]
fps = [ AllChem.GetMorganFingerprintAsBitVect(mol, 2) for mol in mols ]

fps_arr = []
for fp in fps:
    arr = np.zeros((1,), dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fp, arr)
    fps_arr.append(arr)

label_dict = { 'Active':0, 'Inconclusive':1, 'Inactive':1}
labels = [label_dict[mol.GetProp('PUBCHEM_ACTIVITY_OUTCOME')] for mol in mols ]
feats = np.asarray(fps_arr, dtype=np.float32)
labels = np.asarray( labels, dtype=np.int64 ) 

train_feats, val_feats, train_labels, val_labels = train_test_split( feats, labels, test_size=0.2, random_state=100  )

train_obj = gzip.open('ipab_train_bi.pickle.gz', 'wb')
val_obj = gzip.open('ipab_validation_bi.pickle.gz', 'wb' )

cPickle.dump([ train_feats, train_labels ], train_obj, cPickle.HIGHEST_PROTOCOL)
cPickle.dump([ val_feats, val_labels ], val_obj, cPickle.HIGHEST_PROTOCOL)

