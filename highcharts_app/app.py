from flask import Flask, render_template
app = Flask( __name__ )

from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import MolDraw2DSVG
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor, Descriptors, AllChem, DataStructs
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import pickle
#structures from drug bank
drugs = [ mol for mol in Chem.SDMolSupplier( "structures.sdf" ) if mol != None ][:500]
#testset
test = [ mol for mol in Chem.SDMolSupplier( "testset.sdf" ) if mol != None ]


def calc_fp_arr( mols ):
    fplist = []
    for mol in mols:
        arr = np.zeros( (1,) )
        fp = AllChem.GetMorganFingerprintAsBitVect( mol, 2 )
        DataStructs.ConvertToNumpyArray( fp, arr )
        fplist.append( arr )
    return np.asarray( fplist )

def getsvgtext( mol ):
    d2d = rdMolDraw2D.MolDraw2DSVG(200,200)
    d2d.DrawMolecule( mol )
    d2d.FinishDrawing()
    svg = d2d.GetDrawingText()
    return svg.replace( "svg:","" )

drugfparr = calc_fp_arr( drugs )
testfparr = calc_fp_arr( test )

#do PCA
pca = PCA( n_components=2 )
pca.fit( drugfparr )
f = open( 'drugpca.pkl', 'wb' )
pickle.dump( pca, f )
f.close()

drugsX = pca.transform( drugfparr )
data1 = [ { 'x' : drugsX[i][0], 'y':drugsX[i][1], 'svg': getsvgtext( drugs[i] ) } for i in range(len(drugsX)) ]
testX = pca.transform( testfparr )
data2 = [ {  'x': testX[i][0], 'y':testX[i][1], 'svg': getsvgtext( test[i] ) } for i in range(len(testX)) ]


@app.route( '/' )
@app.route( '/chart' )
def chart():
    return render_template( 'chart.html', data1 = data1, data2 = data2 )


if __name__ == '__main__':
    app.debug = True
    app.run(  )
