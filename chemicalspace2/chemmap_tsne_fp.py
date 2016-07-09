import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from sklearn.manifold import TSNE
import pickle, sys
import seaborn as sns
drugs = [ mol for mol in Chem.SDMolSupplier( "gsk3b.sdf" ) if mol != None ]

def calc_fp_arr( mols ):
    fplist = []
    for mol in mols:
        arr = np.zeros( (1,) )
        fp = AllChem.GetMorganFingerprintAsBitVect( mol, 2 )
        DataStructs.ConvertToNumpyArray( fp, arr )
        fplist.append( arr )
    return np.asarray( fplist )

res=calc_fp_arr(drugs)
model = TSNE( n_components=3, init="pca",  n_iter=2000 )
#model = TSNE( n_components=2,random_state=42 ,  n_iter=1000 )
x = model.fit_transform( res )
f = open( 'pca_tsne.pkl', 'wb' )
pickle.dump( model, f )
f.close()



f = open( "pca_tsneres.txt", "w" )
act = []
for i in range( x.shape[0] ):
    line = ""
    line = Chem.MolToSmiles( drugs[i] ) + "," + drugs[i].GetProp( "STANDARD_VALUE" ) + "," + str( x[i][0] )  + "," + str( x[i][1] ) + "\n"
    act.append(  float(drugs[i].GetProp( "STANDARD_VALUE" ))  )
    f.write( line )
f.close()
x = pd.DataFrame( x )
x.columns = [ 'A1', 'A2', 'A3'  ]

g = sns.jointplot( "A1", "A2", data = x )
g.savefig( "TSNE.png" )
