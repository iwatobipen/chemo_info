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
			sssr = Chem.GetSSSR( fragment_obj ) 
			if sssr >= 2:
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
        counter = 1
        
        for mol in sdf:
                try:
                        results = get_bicyclic( mol )
                        for res in results:
                                log.write( "%s\tscaf_%s"%(res,counter) )
                                log.write("\n")
                                counter +=1
                                if counter % 1000 == 0:
                                        print counter
                except:
                        pass

        
        log.close()

        
        
