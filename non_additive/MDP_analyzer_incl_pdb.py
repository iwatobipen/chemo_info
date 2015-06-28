"""
#
# Analyze SAR data for nonadditivity cycles
#
# Copyright 2014 Christian Kramer
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
#     limitations under the License.
#
##############
"""
# This script analyzes matched double pairs and adds pdb information
# Works with python starting from version 2.7
#
# The chembl binding requires chembl to be installed and running in a local mysql database
# RDKit also has to be available and callable
#

import sys, os.path, subprocess, copy, string, math
# sys.path.append('/Users/Kramer/Projects/Matched-Pairs/Code')   # Comment this one out or replace appropriately if you are not Christian Kramer
from MMP_analysis import *
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
        lhs_assay_ids = [cpd["assay_id"] for cpd in meas[moll]]
        rhs_assay_ids = [cpd["assay_id"] for cpd in meas[molr]]
        id_assays = set(lhs_assay_ids).intersection(set(rhs_assay_ids))
        if len(id_assays) == 0: continue            # Cycle if cpds have been measured in different assays
        # Skip pair if the transformation has more than one anchoring point
        # and the topological distance between those changes (no reason to assume additivity)
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
        av_act_r = sum([-math.log10(cpd["value"]/1E9) for cpd in meas[molr] if cpd["assay_id"] in id_assays])/len(id_assays)
        av_act_l = sum([-math.log10(cpd["value"]/1E9) for cpd in meas[moll] if cpd["assay_id"] in id_assays])/len(id_assays)
        if moll in neighs.keys():
            if not molr in [i[1] for i in neighs[moll]]:
                neighs[moll].append((av_act_r,molr,line[4]))
            else:
                idx = [i[1] for i in neighs[moll]].index(molr)
                smirks_len = len(neighs[moll][idx][2])
                if len(line[4]) < smirks_len:
                    neighs[moll][idx] = (av_act_r,molr,line[4])
        else:
            neighs[moll] = [(av_act_r,molr,line[4])]
        if molr in neighs.keys():
            if not moll in [i[1] for i in neighs[molr]]:
                neighs[molr].append((av_act_l,moll,rhs+">>"+lhs))
            else:
                idx = [i[1] for i in neighs[molr]].index(moll)
                smirks_len = len(neighs[molr][idx][2])
                if len(line[4]) < smirks_len:
                    neighs[molr][idx] = (av_act_l,moll,rhs+">>"+lhs)
        else:
            neighs[molr] = [(av_act_l,moll,rhs+">>"+lhs)]

    for moll in neighs.keys():
        nn_mols = sorted(neighs[moll],reverse=True)
        neighs[moll] = [i[1:] for i in nn_mols]

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
    f.write("Circle\tSMIRKS1\tSMIRKS2\tAct1\tAct2\tAct3\tAct4\tDiff1\tDiff2\tDiff3\tDiff4\tAdd_diff\tSimilar_PDBs1\tSimilar_PDBs2\tSimilar_PDBs3\tSimilar_PDBs4\n")

    for i,cpds in enumerate(circles):
        c0_assays = [cpd["assay_id"] for cpd in meas[cpds[0]]]
        c1_assays = [cpd["assay_id"] for cpd in meas[cpds[1]]]
        c2_assays = [cpd["assay_id"] for cpd in meas[cpds[2]]]
        c3_assays = [cpd["assay_id"] for cpd in meas[cpds[3]]]
        id_assays = set(c0_assays).intersection(set(c1_assays)).intersection(set(c2_assays)).intersection(set(c3_assays))
        if len(id_assays) == 0:
            continue
        act0 = sum([-math.log10(cpd["value"]/1E9) for cpd in meas[cpds[0]] if cpd["assay_id"] in id_assays])/len(id_assays)
        act1 = sum([-math.log10(cpd["value"]/1E9) for cpd in meas[cpds[1]] if cpd["assay_id"] in id_assays])/len(id_assays)
        act2 = sum([-math.log10(cpd["value"]/1E9) for cpd in meas[cpds[2]] if cpd["assay_id"] in id_assays])/len(id_assays)
        act3 = sum([-math.log10(cpd["value"]/1E9) for cpd in meas[cpds[3]] if cpd["assay_id"] in id_assays])/len(id_assays)

        acts = [act0,act1,act2,act3]
        mwts = [meas[cpd][0]['mwt'] for cpd in cpds]
        min_idx = mwts.index(min(mwts))
        circles[i] = (circles[i]+circles[i])[min_idx:min_idx+4]
        acts = (acts+acts)[min_idx:min_idx+4]
        act0,act1,act2,act3 = acts

        line = circles[i][0]+"_"+circles[i][1]+"_"+circles[i][2]+"_"+circles[i][3]+"\t"
        idx = [j[0] for j in neighs[circles[i][0]]].index(circles[i][1])
        line = line + neighs[circles[i][0]][idx][1]+"\t"
#        if not "H" in neighs[circles[i][0]][idx][1]: continue
        idx = [j[0] for j in neighs[circles[i][1]]].index(circles[i][2])
        line = line + neighs[circles[i][1]][idx][1]+"\t"
#        if not "H" in neighs[circles[i][1]][idx][1]: continue
        line = line + "{:4.3}".format(act0)+"\t"
        line = line + "{:4.3}".format(act1)+"\t"
        line = line + "{:4.3}".format(act2)+"\t"
        line = line + "{:4.3}".format(act3)+"\t"
        line = line + "{:4.3}".format(act1-act0)+"\t"
        line = line + "{:4.3}".format(act2-act1)+"\t"
        line = line + "{:4.3}".format(act3-act2)+"\t"
        line = line + "{:4.3}".format(act0-act3)+"\t"
        line = line + "{:4.3}".format((act2-act3) - (act1-act0))+"\t"
        line = line + ";".join([k[0]+"_"+"{:3.2}".format(k[1]) for k in pdbs[i][0]])+"\t"
        line = line + ";".join([k[0]+"_"+"{:3.2}".format(k[1]) for k in pdbs[i][1]])+"\t"
        line = line + ";".join([k[0]+"_"+"{:3.2}".format(k[1]) for k in pdbs[i][2]])+"\t"
        line = line + ";".join([k[0]+"_"+"{:3.2}".format(k[1]) for k in pdbs[i][3]])+"\n"
        f.write(line)

    f.close()

    return


def initialize_ChEMBL_PDB_conversion():
    """
    Return ChEMBL to PDB conversion and vice versa via uniprot IDs
    """
    
    f = open("cc-to-pdb.txt",'r')
    a = f.readlines()
    f.close()

    PDB_LIG_ID_to_PDB_PROT = dict((line.split()[0],line.split()[1:]) for line in a)
    PDB_PROT_to_PDB_LIG_ID = defaultdict(list)
    for line in a:
        for pdb_id in line.split()[1:]:
            PDB_PROT_to_PDB_LIG_ID[pdb_id].append(line.split()[0])

    f = open("pdbtosp.txt","r")
    a = f.read()
    f.close()

    a = a.replace("\n                           ","")
    a = a.split("\n")
    a = a[24:-6]

    PDB_RES = dict((string.lower(line.split()[0]),line.split()[2]) for line in a if line.split()[1] == "X-ray")
    PDB_to_UNIPROT = dict((string.lower(line[:4]),[i[:6] for i in line[41:].split("(")]) for line in a if line.split()[1] == "X-ray")

    UNIPROT_to_PDB = defaultdict(list)

    for line in a:
        if not line.split()[1] == "X-ray": continue
        for uniprot_id in [i[:6] for i in line[41:].split("(")]:
            UNIPROT_to_PDB[uniprot_id].append(string.lower(line[:4]))

    # Read PDB Ligand file and generate Fingerprints
    max_heavy = 70

    remover = SaltRemover.SaltRemover()

    pdb_ligands = {}
    f = open("Components-smiles-stereo-oe.smi","r")
    for line in f:
        line = line.split()
        if len(line) >= 2: pdb_ligands[line[1]] = {"smiles":line[0],"fp":None}

    f.close()

    return PDB_LIG_ID_to_PDB_PROT, PDB_PROT_to_PDB_LIG_ID, PDB_RES, PDB_to_UNIPROT, UNIPROT_to_PDB, pdb_ligands


# Get all ChEMBL targets with uniprot code and number of ligands

def get_all_chembl_activities_incl_uniprot_id(activity_type="Ki",user="root",password=""):
    """
        Fetches all entries from CHEMBL17 with specified activity type.
        Uses local mysql installation of ChEMBL. Username and password
        can be supplied using user='' and password=''
        """
    
    up = "-u"+user
    if password != "": up = up + " -p" + password
    
    chembl_string = Template('''
        mysql $up --column-names -e "use chembl_17;
        select activities.assay_id,
        activities.doc_id,
        activities.standard_type,
        activities.standard_relation,
        activities.standard_value,
        activities.standard_units,
        assays.tid,
        target_dictionary.chembl_id,
        component_sequences.accession
        from activities,
        assays,
        target_dictionary,
        component_sequences,
        target_components
        where activities.standard_type='$activity_type' and
        activities.assay_id=assays.assay_id and
        target_components.tid=assays.tid and
        component_sequences.component_id=target_components.component_id and
        assays.tid=target_dictionary.tid;"
        ''')
    
    # Format data
    
    data = subprocess.check_output(chembl_string.substitute(up=up,activity_type=activity_type),shell=True).split("\n")
    data = [line.split() for line in data if len(line) > 0]
    
    return data


def get_target_frequency_count(data,min_activities=3000):
    """
    Extract most frequent targets from extracted chembl data
    """
    header = data[0]
    data = data[1:]
    
    assay_id = header.index("assay_id")
    t_chembl_id = header.index("chembl_id")
    tid = header.index("tid")
    uniprot_id = header.index("accession")
    
    # Remove dummy target
    
    dummy_name = "CHEMBL612545"
    data = [line for line in data if not line[t_chembl_id] == dummy_name]
    
    # Remove empty entries
    
    act_column = header.index("standard_value")
    data = [line for line in data if not line[act_column] == "NULL"]
    
    # Number of entries per assay and unique targets
    
    assay_ids = Counter([i[assay_id] for i in data])
    target_ids = Counter([i[tid] for i in data])
    t_chembl_ids = Counter([i[t_chembl_id] for i in data])
    
    # Most frequent targets
    
    return [(key,value) for key,value in t_chembl_ids.items() if value > min_activities]


def get_pdbs_with_similar_ligands(circles,target,tc,pdbids, PDB_PROT_to_PDB_LIG_ID, pdb_ligands, meas):
    """
    Find similar ligands within pdb files
    """

    remover = SaltRemover.SaltRemover()
    similar_pdbs = []
    for circle in circles:
        circle_pdbs = []
        for cpd in circle:
            t_smi = meas[cpd][0]['smiles']
            t_cpd = Chem.MolFromSmiles(t_smi)
            t_res = remover.StripMol(t_cpd)        # Remove Salts
            t_fp = FingerprintMols.FingerprintMol(t_res)
            pdbs = []
            for pdbid in pdbids:
                for lig in PDB_PROT_to_PDB_LIG_ID[pdbid]:
                    try:
                        if pdb_ligands[lig]["fp"] == "skip": continue
                        if pdb_ligands[lig]["fp"] == None:
                            cpd = Chem.MolFromSmiles(pdb_ligands[lig]["smiles"])
                            res = remover.StripMol(cpd)        # Remove Salts
                            if res.GetNumAtoms() > max_heavy: continue # if the ligand is too large
                            smiles = Chem.MolToSmiles(res)     # Canonicalize smiles
                            fp = FingerprintMols.FingerprintMol(res)
                            pdb_ligands[lig] = ({"smiles":smiles,"fp":fp,"mol":cpd})
                    except:
                        pdb_ligands[lig] = ({"fp":"skip"})
                        continue
                    sim = DataStructs.FingerprintSimilarity(t_fp,pdb_ligands[lig]["fp"])
                    if sim >= tc:
                        pdbs.append((pdbid,sim))
            circle_pdbs.append(pdbs)
        similar_pdbs.append(circle_pdbs)

    return similar_pdbs


def get_circles_with_pdb_info(target,activity_type,mysql_user,mysql_password,min_conf,max_diff,max_heavy,mmpa_home,tc,pdb_ids,pdb_prot_to_lig,pdb_ligands,use_old_results,no_chiral,outfile):
    """
    Extract thermodynamic transformation circles with associated information from ChEMBL and add Info about PDB files with similar ligands.
    """
    data = extract_all_activities_for_target(target,activity_type=activity_type,user=mysql_user,password=mysql_password)
    data = clean_chembl_data(data, min_conf=min_conf)
    meas = build_ligand_dictionary(data, maxdiff=max_diff)
    meas = add_mwt(meas)
        
    if not use_old_results or not os.path.exists(target+"_"+activity_type+"_mmp_raw.csv"):
        write_smiles_id_file(meas,target,activity_type,max_heavy)
        calc_raw_MMPs(target,activity_type,mmpa_home)
        
    mmps = read_raw_mmps(target,activity_type)
    neighs = build_neighbor_dictionary(mmps,meas,no_chiral=no_chiral)
    circles = get_circles(neighs)
    pdbs = get_pdbs_with_similar_ligands(circles,target,tc,pdb_ids,pdb_prot_to_lig,pdb_ligands,meas)
    write_circles_to_output(circles,meas,pdbs,neighs,outfile)


def add_mwt(meas):
    """
    Add MWT information to meas
    """
    from rdkit.Chem import Descriptors
    for entry in meas.keys():
        mol = Chem.MolFromSmiles(meas[entry][0]["smiles"])
        try:
            meas[entry][0]['mwt']=Descriptors.MolWt(mol)
        except:
            del meas[entry]

    return meas



####################
# Defaults

target = ""
mysql_user = "root"
mysql_password = ""
use_old_results = False
activity_type = "Ki"
min_activities = 2000
outfile = ""
max_diff = 2.5
max_heavy = 70
min_conf = 4
smiin = ""
tc = 0.7
no_chiral = False
add_pdb_info = False

Usage = """
    Usage
    python MDP_analyzer.py [-target {id,all}] [-mysql_user id] [-mysql_password id] [-use_old_results]
    [-activity_type {Ki, IC50, both}] [-max_diff #] [-min_Nactivities #] [-out filename] [-tc #]
    [-max_heavy #] [-min_conf #] [-no_chiral] [-out filename] [-add_pdb_info]
    
    Keyword           Default         Meaning
    target            none            CHEMBL Target ID
    mysql_user        root            MYSQL User ID
    mysql_password    none            MYSQL User Password
    use_old_results   False           Use raw MMP results from previous computation
    activity_type     Ki              Select activity type to be used from [Ki, IC50, both]
    max_diff          2.5             Maximum tolerated difference in log(Activity) for multiple
                                      measurements of a given ligand/protein complex.
    max_heavy         70              Maximum number of heavy atoms per ligand
    min_conf          4               Minimum CHEMBL Confidence Score
    min_Nactivities   2000            Minimum Number of activities when calculating cycles for all targets
    smiin             none            SMILES input filename
    tc                0.7             Tanimoto Coefficient for defining similarity to PDB ligands
    no_chiral         False           Skip all transformations that include chiral centers
    add_pdb_info      False           Add Info about PDB entries with similar ligands (works for CHEMBL targets only)

    out               Activity_diffs_target_activity.txt          Output file name
    
    """

# Read input
for i in range(len(sys.argv)):
    if sys.argv[i] == "-target":
        target = sys.argv[i+1]
    elif sys.argv[i] == "-mysql_user":
        mysql_user = sys.argv[i+1]
    elif sys.argv[i] == "-mysql_password":
        mysql_password = sys.argv[i+1]
    elif sys.argv[i] == "-use_old_results":
        use_old_results = True
    elif sys.argv[i] == "-activity_type":
        activity_type = sys.argv[i+1]
    elif sys.argv[i] == "-max_heavy":
        max_heavy = int(sys.argv[i+1])
    elif sys.argv[i] == "-max_diff":
        max_diff = float(sys.argv[i+1])
    elif sys.argv[i] == "-out":
        outfile = sys.argv[i+1]
    elif sys.argv[i] == "-smiin":
        smiin = sys.argv[i+1]
    elif sys.argv[i] == "-min_conf":
        min_conf = int(sys.argv[i+1])
    elif sys.argv[i] == "-min_Nactivities":
        min_activities = int(sys.argv[i+1])
    elif sys.argv[i] == "-tc":
        tc = float(sys.argv[i+1])
    elif sys.argv[i] == "-no_chiral":
        no_chiral = True
    elif sys.argv[i] == "-add_pdb_info":
        add_pdb_info = True
    elif sys.argv[i] == "-h":
        print Usage
        exit(0)


mmpa_home = subprocess.check_output("echo $RDBASE",shell=True).split()[0]+"/Contrib/mmpa/"

if len(sys.argv) < 2:
    print Usage
    exit(0)

if not target[:6] == "CHEMBL" or target == "all" or smiin != "":
    print "Target not recognized. Either supply SMILES input file or target named 'CHEMBL...' or 'all'."
    print "Exiting."
    print Usage
    exit(0)

if activity_type not in ["Ki","IC50","both"] and smiin == "":
    print "Allowed activity types for analyzing CHEMBL data: Ki, IC50, both"
    print "Exiting"
    exit(0)

if smiin != "":
    if add_pdb_info:
        print "Adding PDB Info to ligands from SMILES input files currently not supported."
    meas = build_ligand_dictionary_from_infile(smiin)
    meas = add_mwt(meas)
    if use_old_results and not os.path.exists(target+"_"+activity_type+"_mmp_raw.csv"): use_old_results = False
    
    if not use_old_results:
        write_smiles_id_file(meas,target,activity_type,max_heavy)
        calc_raw_MMPs(target,activity_type,mmpa_home)
    
    mmps = read_raw_mmps(target,activity_type)
    neighs = build_neighbor_dictionary(mmps,meas,no_chiral=no_chiral)
    circles = get_circles(neighs)
    write_circles_to_output(circles,meas,outfile)
    exit(0)

PDB_LIG_ID_to_PDB_PROT, PDB_PROT_to_PDB_LIG_ID, PDB_RES, PDB_to_UNIPROT, UNIPROT_to_PDB, pdb_ligands = initialize_ChEMBL_PDB_conversion()

if activity_type == "both":
    ki_data = get_all_chembl_activities_incl_uniprot_id(activity_type="Ki")
    ic50_data = get_all_chembl_activities_incl_uniprot_id(activity_type="IC50")
    data = ki_data + ic50_data[1:]
    ki_data = ""
    ic50_data = ""
    CHEMBL_ID_to_UNIPROT = dict((line[-2],line[-1]) for line in data)
    UNIPROT_to_CHEMBL_ID = dict((line[-1],line[-2]) for line in data)
else:
    data = get_all_chembl_activities_incl_uniprot_id(activity_type=activity_type)
    CHEMBL_ID_to_UNIPROT = dict((line[-2],line[-1]) for line in data)
    UNIPROT_to_CHEMBL_ID = dict((line[-1],line[-2]) for line in data)

if target == "all" and activity_type == "both":
    most_frequent_targets = get_target_frequency_count(data,min_activities)
    for entry in most_frequent_targets:
        target = entry[0]
        upid = CHEMBL_ID_to_UNIPROT[target]
        pdbids = UNIPROT_to_PDB[upid]
        if len(pdbids) == 0: continue
        ligs = [PDB_PROT_to_PDB_LIG_ID[pdbid] for pdbid in pdbids]
        if len([i for i in j for j in ligs]) == 0: continue
        # Ki data first
        activity_type = "Ki"
        outfile = "Additivity_diffs_"+target+"_Ki.txt"
        get_circles_with_pdb_info(target,activity_type,mysql_user,mysql_password,min_conf,max_diff,max_heavy,mmpa_home,tc,UNIPROT_to_PDB[CHEMBL_ID_to_UNIPROT[target]],PDB_PROT_to_PDB_LIG_ID,pdb_ligands,use_old_results,no_chiral,outfile)
        # IC50 data next
        activity_type = "IC50"
        outfile = "Additivity_diffs_"+target+"_IC50.txt"
        get_circles_with_pdb_info(target,activity_type,mysql_user,mysql_password,min_conf,max_diff,max_heavy,mmpa_home,tc,UNIPROT_to_PDB[CHEMBL_ID_to_UNIPROT[target]],PDB_PROT_to_PDB_LIG_ID,pdb_ligands,use_old_results,no_chiral,outfile)
        exit(0)
elif target == "all":
    most_frequent_targets = get_target_frequency_count(data,min_activities)
    for entry in most_frequent_targets:
        target = entry[0]
        upid = CHEMBL_ID_to_UNIPROT[target]
        pdbids = UNIPROT_to_PDB[upid]
        if len(pdbids) == 0: continue
        ligs = [PDB_PROT_to_PDB_LIG_ID[pdbid] for pdbid in pdbids]
        if len([i for i in j for j in ligs]) == 0: continue
        outfile = "Additivity_diffs_"+target+"_"+activity_type+".txt"
        get_circles_with_pdb_info(target,activity_type,mysql_user,mysql_password,min_conf,max_diff,max_heavy,mmpa_home,tc,UNIPROT_to_PDB[CHEMBL_ID_to_UNIPROT[target]],PDB_PROT_to_PDB_LIG_ID,pdb_ligands,use_old_results,no_chiral,outfile)
        exit(0)
elif activity_type == "both":
    # Ki data first
    activity_type = "Ki"
    outfile = "Additivity_diffs_"+target+"_Ki.txt"
    get_circles_with_pdb_info(target,activity_type,mysql_user,mysql_password,min_conf,max_diff,max_heavy,mmpa_home,tc,UNIPROT_to_PDB[CHEMBL_ID_to_UNIPROT[target]],PDB_PROT_to_PDB_LIG_ID,pdb_ligands,use_old_results,no_chiral,outfile)
    # IC50 data next
    activity_type = "IC50"
    outfile = "Additivity_diffs_"+target+"_IC50.txt"
    get_circles_with_pdb_info(target,activity_type,mysql_user,mysql_password,min_conf,max_diff,max_heavy,mmpa_home,tc,UNIPROT_to_PDB[CHEMBL_ID_to_UNIPROT[target]],PDB_PROT_to_PDB_LIG_ID,pdb_ligands,use_old_results,no_chiral,outfile)
    exit(0)
else:
    if outfile == "": outfile = "Additivity_diffs_"+target+"_"+activity_type+".txt"
    get_circles_with_pdb_info(target,activity_type,mysql_user,mysql_password,min_conf,max_diff,max_heavy,mmpa_home,tc,UNIPROT_to_PDB[CHEMBL_ID_to_UNIPROT[target]],PDB_PROT_to_PDB_LIG_ID,pdb_ligands,use_old_results,no_chiral,outfile)
    exit(0)

