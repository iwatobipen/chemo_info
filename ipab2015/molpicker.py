from rdkit import Chem
cids = open('selectid.txt', 'r')
cids = [ line.rsplit()[0] for line in cids  ]

sdf = Chem.ForwardSDMolSupplier('ContestLibrary.sdf')
#print [ i for i in mol.GetPropNames()]
#print mol.GetProp('idnumber')

outf = open('selectedmols.smi', 'w')
counter = 0
for mol in sdf:
    try:
        cid = mol.GetProp('idnumber')
	if cid in cids:
		smi = Chem.MolToSmiles(mol)
		outf.write(smi+'\t'+cid+'\n')
	counter +=1
	if counter % 50000 ==0:
		print counter
    except:
        pass
outf.close()

