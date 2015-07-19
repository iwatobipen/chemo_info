import gzip
import time
import cPickle
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import numpy as np
f = open('timelog.txt', 'w')
logger = open('indx.txt', 'w')
sdf = Chem.ForwardSDMolSupplier('ContestLibrary.sdf')
fps_arr = [] 
count = 0
file_idx = 1
lap = 0
start = time.time()
for mol in sdf:
    if mol is not None:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
        arr = np.zeros((1,), dtype=np.float32)
        DataStructs.ConvertToNumpyArray( fp, arr )
        fps_arr.append( arr )
        idn = mol.GetProp('idnumber')
        logger.write('%s\n'%idn)
        count += 1
        if count % 10000 == 0:
            old_lap = lap
            lap = time.time()
            total_delta = lap-start
            lap_delta = lap-old_lap
            f.write('%s mols\t%s\t%s\n'%(count, total_delta, lap_delta))
        if count % 100000 == 0:
            print count
            fname = './contest_feat/contestlib_%s.pickle.gz'%file_idx
            file_obj = gzip.open(fname, 'wb')
            cPickle.dump( [np.asarray(fps_arr), np.zeros(len(fps_arr),dtype=np.int64)], file_obj, cPickle.HIGHEST_PROTOCOL )
            file_idx += 1
            file_obj.close()
            fps_arr = []
f.close()
file_idx += 1
file_obj = gzip.open('./contest_feat/contestlib_%s.pickle.gz'%file_idx, 'wb')
cPickle.dump( [np.asarray(fps_arr), np.zeros((len(fps_arr),))], file_obj, cPickle.HIGHEST_PROTOCOL )
file_obj.close()
logger.close()

