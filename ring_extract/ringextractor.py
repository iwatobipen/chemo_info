import sys
from rdkit import Chem
from rdkit.Chem import Recap


def get_bicyclic( mol ):
        bicyclic = []
	leaves = Recap.RecapDecompose( mol ).GetLeaves()
	if len( leaves ) != 0:
		fragments = leaves.keys()
		for fragment in fragments:
			fragment_obj = Chem.MolFromSmiles( fragment )
			ring_info = fragment_obj.GetRingInfo()
			if ring_info.NumRings() >= 2:
				scaffold = Chem.MurckoDecompose( fragment_obj )
				bicyclic.append( Chem.MolToSmiles( scaffold ) )
	else:
                sssr = Chem.GetSSSR(mol)
                if sssr >= 2:
                        scaffold = Chem.MurckoDecompose( mol )
                        bicyclic.append( Chem.MolToSmiles( scaffold ) )
        return bicyclic


if __name__ == '__main__':
        sdf = Chem.SDMolSupplier( sys.argv[1] )
        log = open( "rings.txt", "w" )
        scaf_set = set()
        for mol in sdf:
                try:
                        results = get_bicyclic( mol )
                        for res in results:
                                scaf_set.add( res )                      
                except:
                        pass
        sdf = []
        
        for idx,scaf in enumerate(scaf_set):
                log.write("%s\tscaf_%s"%(scaf,idx))
                log.write("\n")
        log.close()

        
        
