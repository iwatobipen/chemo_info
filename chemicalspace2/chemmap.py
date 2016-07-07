import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from sklearn.decomposition import PCA
import pickle, sys

drugs = [ mol for mol in Chem.SDMolSupplier( "drugs.sdf" ) if mol != None ]

def calc_fp_arr( mols ):
    fplist = []
    for mol in mols:
        arr = np.zeros( (1,) )
        fp = AllChem.GetMorganFingerprintAsBitVect( mol, 2 )
        DataStructs.ConvertToNumpyArray( fp, arr )
        fplist.append( arr )
    return np.asarray( fplist )

res=calc_fp_arr(drugs)
pca = PCA( n_components = 3 )
pca.fit( res )
f = open( 'pca.pkl', 'wb' )
pickle.dump( pca, f )
f.close()


x = pca.transform( res )
f = open( "pcares.txt", "w" )

for i in range( x.shape[0] ):
    line = ""
    line = Chem.MolToSmiles( drugs[i] ) + "," + str( x[i][0] )  + "," + str( x[i][1] ) + "\n"
    f.write( line )
f.close()
