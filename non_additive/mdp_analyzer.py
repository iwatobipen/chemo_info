import sys, os.path, subprocess, copy, string, math
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import SaltRemover
from rdkit.Chem.Fingerprints import FingerprintMols

from string import Template
from collections import Counter
from collections import defaultdict

mmps = open("sample_mmps_default.csv","r")
mmps = [line.rstrip() for line in mmps]
def build_neighbor_dictionary( mmps, no_chiral=False ):
    """
    Generate neighbor dictionary from mmps read
    """
    neighs = {}
    for line in mmps:
        line = line.split( "," )
        if no_chiral:
            if "@" in line[ 4 ]: continue
        lhs = line[4].split(">>")[0]
        rhs = line[4].split(">>")[1]
        moll = line[0]
        molr = line[1]
        # skip pair if the transformation has more than one anchoring point
        # and the topological distance between those changes
        if "[*:2]" in lhs:
            a = Chem.MolFromSmarts(lhs)
            b = Chem.MolFromSmarts(rhs)
            a_idx1 = [ atom.GetSmarts() for atom in a.GetAtoms() ].index("[*:1]")
            a_idx2 = [ atom.GetSmarts() for atom in a.GetAtoms() ].index("[*:2]")
            b_idx1 = [ atom.GetSmarts() for atom in b.GetAtoms() ].index("[*:1]")
            b_idx2 = [ atom.GetSmarts() for atom in b.GetAtoms() ].index("[*:2]")
            if not Chem.GetDistanceMatrix(a)[a_idx1,a_idx2] == Chem.GetDistanceMatrix(b)[b_idx1,b_idx2]: continue
            if "[*:3]" in lhs:
                a_idx3 = [ atom.GetSmarts() for atom in a.GetAtoms() ].index("[*:3]")
                b_idx3 = [ atom.GetSmarts() for atom in b.GetAtoms() ].index("[*:3]")
                if not Chem.GetDistanceMatrix(a)[a_idx1,a_idx3] == Chem.GetDistanceMatrix(b)[b_idx1,b_idx3]: continue
                if not Chem.GetDistanceMatrix(a)[a_idx1,a_idx3] == Chem.GetDistanceMatrix(b)[b_idx1,b_idx3]: continue
        if moll in neighs.keys():
            if not molr in [i[1] for i in neighs[moll]]:
                neighs[moll].append((molr,line[4]))
            else:
                idx = [i[1] for i in neighs[moll]].index(molr)
                smirks_len = len( neighs[moll][idx][2] )
                if len(line[4]) < smirks_len:
                    nighs[moll][idx] = (molr, line[4])
        else:
            neighs[ moll ] = [(molr, line[4])]
        if molr in neighs.keys():
            if not molr in [i[1] for i in neighs[molr]]:
                neighs[molr].append((moll, rhs+">>"+lhs))
            else:
                idx = [i[1] for i in neighs[molr]].index(moll)
                smirks_len = len(neighs[molr][idx][2])
                if len(line[4]) < smirks_len:
                    neighs[molr][idx] = (moll, rhs+">>"+lhs)
        else:
            neighs[molr] = [(moll, rhs+">>"+lhs)]
    for moll in neighs.keys():
        nn_mols = sorted(neighs[moll], reverse=True)
        #neighs[moll] = [i[1:] for i in nn_mols]
        neighs[moll] = [i for i in nn_mols]
    return neighs

def get_circles(neighs):
    """
    Assemble circle information
    """
    neighs1 = copy.deepcopy(neighs)
   
    circles = []
    for k1 in neighs1.keys():
        neighs2 = copy.deepcopy(neighs1)
        for k2 in neighs1[k1]:
            if k2[0] in neighs1.keys():
                tr1_back = k2[1].split(">>")[1]+">>"+k2[1].split(">>")[0]
                for k3 in neighs1[k2[0]]:
                    if k3[0] == k1: continue
                    if k3[0] in neighs1.keys():
                        if k1 in [i[0] for i in neighs[k3[0]]]: continue
                        tr2_back = k3[1].split(">>")[1]+">>"+k3[1].split(">>")[0]
                        for k4 in neighs1[k3[0]]:
                            if not k4[1] == tr1_back: continue
                            if k4[0] == k2[0]: continue
                            if k4[0] in neighs2.keys():
                                for k5 in neigh1[k4[0]]:
                                    if not k5[1] == tr2_back: continue
                                    if k5[0] == k1: circles.append((k1,k2[0],k3[0],k4[0]))
                    del neighs2[k2[0]]
        del neighs1[k1]
    return circles


n = build_neighbor_dictionary(mmps)
c = get_circles(n)
print c
