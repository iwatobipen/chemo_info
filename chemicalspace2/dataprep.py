import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

df = pd.read_table( "drugs.txt", header=0 )
df = df[ df.smiles != None ]
rows = df.shape[ 0 ]
mols = []

sd = Chem.SDWriter( "drugs.sdf" )
for i in range( rows ):
    try:
        mol = Chem.MolFromSmiles( df.smiles[ i ] )
        mol.SetProp( "cns_drug", str( df.cns_drug[ i ] ))
        mol.SetProp( "generic_name", df.generic_name[ i ] )
        sd.write( mol )
    except:
        pass
sd.close()
