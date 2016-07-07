import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from sklearn.manifold import TSNE
import pickle, sys

drugs = [ mol for mol in Chem.SDMolSupplier( "st178.sdf" ) if mol != None ]

def calc_fp_arr( mols ):
    fplist = []
    for mol in mols:
        arr = np.zeros( (1,) )
        fp = AllChem.GetMorganFingerprintAsBitVect( mol, 2 )
        DataStructs.ConvertToNumpyArray( fp, arr )
        fplist.append( arr )
    return np.asarray( fplist )

res=calc_fp_arr(drugs)
#model = TSNE( n_components=2,init="pca",  n_iter=1000 )
model = TSNE( n_components=2,random_state=42 ,  n_iter=1000 )
x = model.fit_transform( res )
f = open( 'rand_tsne.pkl', 'wb' )
pickle.dump( model, f )
f.close()


f = open( "rand_tsneres.txt", "w" )

for i in range( x.shape[0] ):
    line = ""
    line = Chem.MolToSmiles( drugs[i] ) + "," + str( x[i][0] )  + "," + str( x[i][1] ) + "\n"
    f.write( line )
f.close()
