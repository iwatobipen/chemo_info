import sys, os.path, subprocess, copy, string, math
# sys.path.append('/Users/Kramer/Projects/Matched-Pairs/Code')   # Comment this one out or replace appropriately if you are not Christian Kramer
from string import Template
from collections import Counter
from collections import defaultdict
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import SaltRemover
from rdkit.Chem.Fingerprints import FingerprintMols



def build_neighbor_dictionary(mmps,meas,no_chiral=False):
    """
    Generate neighbor dictionary from mmps read
    """

    neighs = {}
    for line in mmps:
        line = line.split(",")
        if no_chiral:
            if "@" in line[4]: continue
        lhs = line[4].split(">>")[0]
        rhs = line[4].split(">>")[1]
#        if not "[H]" in line[4]: continue    # Filter out pairs where the transformation exchanges too large fragments
        moll = line[2]
        molr = line[3]

        if "[*:2]" in lhs:
            a = Chem.MolFromSmarts(lhs)
            b = Chem.MolFromSmarts(rhs)
            a_idx1 = [atom.GetSmarts() for atom in a.GetAtoms()].index("[*:1]")
            a_idx2 = [atom.GetSmarts() for atom in a.GetAtoms()].index("[*:2]")
            b_idx1 = [atom.GetSmarts() for atom in b.GetAtoms()].index("[*:1]")
            b_idx2 = [atom.GetSmarts() for atom in b.GetAtoms()].index("[*:2]")
            if not Chem.GetDistanceMatrix(a)[a_idx1,a_idx2] == Chem.GetDistanceMatrix(b)[b_idx1,b_idx2]: continue
            if "[*:3]" in lhs:
               a_idx3 = [atom.GetSmarts() for atom in a.GetAtoms()].index("[*:3]")
               b_idx3 = [atom.GetSmarts() for atom in b.GetAtoms()].index("[*:3]")
               if not Chem.GetDistanceMatrix(a)[a_idx1,a_idx3] == Chem.GetDistanceMatrix(b)[b_idx1,b_idx3]: continue
               if not Chem.GetDistanceMatrix(a)[a_idx2,a_idx3] == Chem.GetDistanceMatrix(b)[b_idx2,b_idx3]: continue
        # Add to neighbor dictionary
        if moll in neighs.keys():
            if not molr in [i[0] for i in neighs[moll]]:
                neighs[moll].append((molr,line[4]))
            else:
                idx = [i[0] for i in neighs[moll]].index(molr)
                smirks_len = len(neighs[moll][idx][1])
                if len(line[4]) < smirks_len:
                    neighs[moll][idx] = (molr,line[4])
        else:
            neighs[moll] = [(molr,line[4])]
        if molr in neighs.keys():
            if not moll in [i[0] for i in neighs[molr]]:
                neighs[molr].append((moll,rhs+">>"+lhs))
            else:
                idx = [i[0] for i in neighs[molr]].index(moll)
                smirks_len = len(neighs[molr][idx][1])
                if len(line[4]) < smirks_len:
                    neighs[molr][idx] = (moll,rhs+">>"+lhs)
        else:
            neighs[molr] = [(moll,rhs+">>"+lhs)]

    for moll in neighs.keys():
        nn_mols = sorted(neighs[moll],reverse=True)
        neighs[moll] = [i[:] for i in nn_mols]

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
                                for k5 in neighs1[k4[0]]:
                                    if not k5[1] == tr2_back: continue
                                    if k5[0] == k1: circles.append((k1,k2[0],k3[0],k4[0]))
                del neighs2[k2[0]]
        del neighs1[k1]

    return circles


def write_circles_to_output(circles,meas,pdbs,neighs,outfile):
    """
    Get Diffs and write to output
    """

    f = open(outfile,'w')
    #f.write("Circle\tSMIRKS1\tSMIRKS2\tAct1\tAct2\tAct3\tAct4\tDiff1\tDiff2\tDiff3\tDiff4\tAdd_diff\tSimilar_PDBs1\tSimilar_PDBs2\tSimilar_PDBs3\tSimilar_PDBs4\n")
    f.write("Circle\tSMIRKS1\tSMIRKS2\n")


    for i,cpds in enumerate(circles):

        mwts = [meas[cpd][0]['mwt'] for cpd in cpds]
        min_idx = mwts.index(min(mwts))
        circles[i] = (circles[i]+circles[i])[min_idx:min_idx+4]

        line = circles[i][0]+"_"+circles[i][1]+"_"+circles[i][2]+"_"+circles[i][3]+"\t"
        idx = [j[0] for j in neighs[circles[i][0]]].index(circles[i][1])
        line = line + neighs[circles[i][0]][idx][1]+"\t"
#        if not "H" in neighs[circles[i][0]][idx][1]: continue
        idx = [j[0] for j in neighs[circles[i][1]]].index(circles[i][2])
        line = line + neighs[circles[i][1]][idx][1]+"\n"
        f.write(line)

    f.close()

    return

