from rdkit import Chem
from rdkit.Chem import ChemicalFeatures # 2D PharmacoPhore FingerPrints
from rdkit.Chem import Pharm2D
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.Pharm2D import SigFactory
from rdkit.Chem.Pharm2D import DefaultSigFactory
from rdkit import RDConfig
from rdkit import DataStructs
import os, sys

dataDir  = os.path.join(RDConfig.RDCodeDir,'Chem/Pharm2D/test_data')
featFact = ChemicalFeatures.BuildFeatureFactory(os.path.join(dataDir,'BaseFeatures.fdef'))
print "########"
sigFact = SigFactory.SigFactory( featFact )
bins = tuple([(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,100)])
sigFact.SetBins(bins)
sigFact.Init()


def calc_p4_fp( mol ):
    fp = Generate.Gen2DFingerprint( mol, sigFact )
    return fp

def calc_tanimoto ( fp1, fp2 ):
    tc = DataStructs.TanimotoSimilarity( fp1,fp2 )
    return tc

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "You need 2 arguments, input.sdf and output.txt"
        

    mols = [ mol for mol in Chem.SDMolSupplier(sys.argv[1]) ]

    print str(len(mols)) + "read."

    output = open( sys.argv[2], "w" )
    output.write( "mol1\tmol2\tp4core_tc" )
    output.write( "\n" )

    fps = {}
    error = {}
    for i,mol in enumerate(mols):
        try:
            fp = calc_p4_fp( mol )
            fps[i] = fp
        except:
            error[i] = Chem.MolToSmiles( mol )

    
    for i in range( len(fps) ):
        for j in range( i ):
            mol1_smi = Chem.MolToSmiles( mols[fps.keys()[i]] )
            mol2_smi = Chem.MolToSmiles( mols[fps.keys()[j]] )
            p4_tanimoto = calc_tanimoto( fps[fps.keys()[i]], fps[fps.keys()[j]] )
            output.write( "%s\t%s\t%s"%( mol1_smi, mol2_smi, p4_tanimoto ) )
            output.write( "\n" )
    output.close()

    error_log = open("errormol.txt","w")
    for i,smi in error.items():
        error_log.write("%s\t%s\n"%(i, smi))
    error_log.close()
    
