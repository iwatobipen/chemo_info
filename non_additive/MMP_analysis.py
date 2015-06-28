"""
#
# Analyze ChEMBL data content which is suitable for MMP analysis
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
# Get chembl data
"""
#pylint: disable=W0311
import math, subprocess
import os.path
from string import Template
from scipy import stats
from collections import Counter


def get_all_chembl_activities(activity_type="Ki",user="root",password=""):
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
            target_dictionary.chembl_id
       from activities,
            assays,
            target_dictionary
      where activities.standard_type='$activity_type' and
            activities.assay_id=assays.assay_id and
            assays.tid=target_dictionary.tid;"
    ''')

# Format data

    data = subprocess.check_output(chembl_string.substitute(up=up,activity_type=activity_type),shell=True).split("\n")
    data = [line.split() for line in data if len(line) > 0]

    return data

def get_most_frequent_targets(data,min_activities=3000):
    """
Extract most frequent targets from extracted chembl data
    """
    header = data[0]
    data = data[1:]

    assay_id = header.index("assay_id")
    t_chembl_id = header.index("chembl_id")
    tid = header.index("tid")

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

    return [key for key in t_chembl_ids.keys() if t_chembl_ids[key]> min_activities]


def extract_all_activities_for_target(target,activity_type="Ki",user="root",password=""):
    """
Extract all activity data for specific CHEMBL target ID
    """
    if target == "":
        print "A target name has to be supplied. Exiting"
        exit(0)
    
    up = "-u"+user
    if password != "": up = up + " -p" + password

    mysql_string = Template('''
    mysql $up --column-names -e "use chembl_17;
    select assays.chembl_id,
           activities.standard_type,
           activities.standard_relation,
           activities.standard_value,
           activities.standard_units,
           target_dictionary.chembl_id,
           compound_structures.canonical_smiles,
           molecule_dictionary.chembl_id,
           docs.chembl_id,
           assays.confidence_score
    from   activities,
           assays,
           compound_structures,
           molecule_dictionary,
           target_dictionary,
           docs
    where  activities.standard_type='$activity_type' and
           target_dictionary.chembl_id='$target' and
           activities.assay_id=assays.assay_id and
           activities.molregno = compound_structures.molregno and
           activities.molregno = molecule_dictionary.molregno and
           activities.doc_id=docs.doc_id and
           assays.tid=target_dictionary.tid;"
    ''')

    print "Extracting data from ChEMBL"
    sdata = subprocess.check_output(mysql_string.substitute(target=target,up=up,activity_type=activity_type),shell=True).split("\n")
    sdata = [line.split() for line in sdata if len(line) > 0]

    return sdata


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def clean_chembl_data(data, min_conf=4, l_act_t = 0.0001, u_act_t = 1000000):
    """
    Clean up data downloaded directly from ChEMBL
    """
    header = data[0]
    
    act_col = header.index("standard_value")
    rel_col = header.index("standard_relation")
    smi_col = header.index("canonical_smiles")
    unit_col = header.index("standard_units")
    conf_col = header.index("confidence_score")
    
    # Remove empty entries
    data = [line for line in data if not line[act_col] == "NULL"]

    # Remove entries that don't have a number in the activity column
    data = [line for line in data if is_number(line[act_col])]

    # Remove entries with nonsense values
    data = [line for line in data if float(line[act_col]) < u_act_t]
    data = [line for line in data if float(line[act_col]) > l_act_t]
    
    # Remove entries with qualified activities
    data = [line for line in data if line[rel_col] == "="]
    
    # Remove entries without smiles
    data = [line for line in data if not "CHEMBL" in line[smi_col]]
    
    # Remove entries that do not have the standard unit "nM"
    data = [line for line in data if line[unit_col] == "nM"]

    # Remove entries with too low confidence score
    data = [line for line in data if is_number(line[conf_col])]
    data = [line for line in data if int(line[conf_col]) >= min_conf]

    data = [header] + data
    return data


def build_ligand_dictionary(data, maxdiff=2.5):
    """
    Return dictionary with measured values
    and other info per ligand.
    maxdiff is the maximum difference between measured value per ligand.
    If the difference is larger, the whole ligand entry is deleted.
    """

    print "Building ligand dictionary"

    header = data[0]
    header[header.index('chembl_id')] = "a_chembl_id"
    header[header.index('chembl_id')] = "t_chembl_id"
    header[header.index('chembl_id')] = "m_chembl_id"
    header[header.index('chembl_id')] = "d_chembl_id"

    assay_id = header.index("a_chembl_id")
    cpd_id = header.index("m_chembl_id")
    act_col = header.index("standard_value")
    rel_col = header.index("standard_relation")
    smi_col = header.index("canonical_smiles")
    unit_col = header.index("standard_units")
    doc_id = header.index("d_chembl_id")
    conf_col = header.index("confidence_score")
    data = data [1:]

# Get unique assays and molecules
    assay_ids = Counter([i[assay_id] for i in data])
    cpd_ids = Counter([i[cpd_id] for i in data])

# Assemble dictionary per compound
    meas = dict()
    for line in data:
        if line[cpd_id] in meas.keys():
            meas[line[cpd_id]].append(dict(assay_id=line[assay_id],value=float(line[act_col]),unit=line[unit_col],smiles=line[smi_col].replace("\\\\","\\"),doc=line[doc_id]))
        else:
            meas[line[cpd_id]]= [dict(assay_id=line[assay_id],value=float(line[act_col]),unit=line[unit_col],smiles=line[smi_col].replace("\\\\","\\"),doc=line[doc_id])]

# Clean dataset from multiple entries
# 1st: remove identical double measurements (keep the entries which originate from assays with more compounds)
    for mol in meas.keys():
        if len(meas[mol]) == 1: continue
        values = [entry["value"] for entry in meas[mol]]
        if len(set(values)) == len(values): continue
        multiples = {}
        for entry in meas[mol]:
            if entry["value"] in multiples.keys():
                multiples[entry["value"]].append(entry)
            else:
                multiples[entry["value"]] = [entry]
        new_list = []
        for entries in multiples.values():
            if len(entries) == 1:
                new_list.append(entries[0])
            else:
                same_assay_counts = []
                for entry in entries:
                    same_assay_counts.append(assay_ids[entry["assay_id"]])
                new_list.append(entries[same_assay_counts.index(max(same_assay_counts))])
        meas[mol] = new_list


# 2nd for the remaining compounds: remove if max and min are too far apart
    for mol in meas.keys():
        if len(meas[mol]) == 1: continue
        values = [entry["value"] for entry in meas[mol]]
        if abs(-math.log10(max(values)/1E9) + math.log10(min(values)/1E9)) > maxdiff:
            del meas[mol]

# 3rd use average if there are several values from identical assays and identical publications available
    for mol in meas.keys():
        if len(meas[mol]) == 1: continue
        assay_ids = [entry["assay_id"] for entry in meas[mol]]
        if len(assay_ids) == len(set(assay_ids)): continue
        assay_ids = list(set(assay_ids))
        n_assay_ids = Counter([entry["assay_id"] for entry in meas[mol]])
        av_act = [0 for i in assay_ids]
        for entry in meas[mol]:
            idx = assay_ids.index(entry["assay_id"])
            av_act[idx] += math.log10(entry["value"]/1E9)
        for i in xrange(len(av_act)): av_act[i] = av_act[i]/n_assay_ids[assay_ids[i]]
        new_list = ["" for i in assay_ids]
        for entry in meas[mol]:
            idx = assay_ids.index(entry["assay_id"])
            new_list[idx] = entry
            new_list[idx]["value"] = 1E9*10**av_act[idx]
        meas[mol] = new_list

    return meas


def write_smiles_id_file(meas,target,activity_type,max_heavy):
    """
Format dataset for MMP analysis and write to temp file
    """

    from rdkit import Chem
    from rdkit.Chem import SaltRemover
    remover = SaltRemover.SaltRemover()

    smifile = target+"_"+activity_type+"_ligands.smi"
    f = open(smifile,'w')

    error_files = target+"_"+activity_type+"_problem_smiles.smi"
    g = open(error_files,'w')

    for mol in meas.keys():
        try:
            cpd = Chem.MolFromSmiles(meas[mol][0]['smiles']) 
            res = remover.StripMol(cpd)        # Remove Salts
            if res.GetNumAtoms() > max_heavy: continue
            smiles = Chem.MolToSmiles(res,True)     # Canonicalize smiles
            if "." in smiles:
                print "Found unknown salt in ",mol,": ",smiles
                print "This compound will be ignored in all further calculations."
                continue

            f.write(smiles+" "+mol+'\n')
        except:
            g.write(meas[mol][0]['smiles']+" "+mol+"\n")

    f.close()
    g.close()

    return


def calc_raw_MMPs(target,activity_type,mmpa_home="/usr/local/share/RDKit/Contrib/mmpa/"):
    """
Generate MMP Indexing and Matching using Jameed Hussain's Code
    """
    smifile = target+"_"+activity_type+"_ligands.smi"
    fragfile = target+"_"+activity_type+"_fragments.txt"
    mmp_outfile = target+"_"+activity_type+"_mmp_raw.csv"

    print "Generating MMP Fragments for "+target
    sp_call = "python "+mmpa_home+"rfrag.py < "+smifile+" > "+fragfile
    subprocess.call(sp_call,shell=True)        # Indexing

    print "Indexing MMP Fragments for "+target+"\n"
    sp_call = "python "+mmpa_home+"indexing.py <"+fragfile+" >"+mmp_outfile
    subprocess.call(sp_call,shell=True)        # Fragment_matching


def read_raw_mmps(target,activity_type):
    """
Read raw precalculated MMPs
    """
    mmp_outfile = target+"_"+activity_type+"_mmp_raw.csv"
    f = open(mmp_outfile,'r')
    mmps = f.readlines()
    f.close()

    return mmps


def workup_raw_transformations(mmps,meas,keep_id_diffs=False):
    """
Clean and unite raw transformations found using Jameed Hussain's MMP code.
    """
# Unite ligand pairs per transformation
# Transformation SMIRKS are in 5th column, pair IDs in 3rd and 4th column

    print "Working up raw transformations"
    transforms = {}
    for line in mmps:
        line = line.split(',')
        transf = line[4].split(">>")
        rev = False
        if len(transf[0]) > len(transf[1]): rev = True
        if transf[1].count("H") > transf[0].count("H"): rev = True
        if not rev:
            transf = line[4]
            pair = [line[2],line[3]]
        else:
            transf = line[4].split(">>")[1]+">>"+line[4].split(">>")[0]
            pair = [line[3],line[2]]
        if transf in transforms.keys():
            if not pair in transforms[transf]:
                transforms[transf].append(pair)
        else:
            transforms[transf] = [pair]

# Add Activity and Diff Information for each pair
# Try to match assay IDs and join from the same assay. If thats not possible, use average activities
    for tf in transforms.keys():
        activities1 = []
        activities2 = []
        identical_assay= []
        deltas = []
        for pair in transforms[tf]:
            assay_ids1 = [entry["assay_id"] for entry in meas[pair[0]]]
            assay_ids2 = [entry["assay_id"] for entry in meas[pair[1]]]
            overlap = list(set(assay_ids1).intersection(set(assay_ids2)))
            if len(overlap) == 0:    # Maybe remove pairs if difference between max and min are too large here
                same_assay = False
                activity1 = sum([-math.log10(float(entry["value"])/1E9) for entry in meas[pair[0]]])/len([entry["value"] for entry in meas[pair[0]]])
                activity2 = sum([-math.log10(float(entry["value"])/1E9) for entry in meas[pair[1]]])/len([entry["value"] for entry in meas[pair[1]]])
                activities1.append(activity1)
                activities2.append(activity2)
                deltas.append(activity2-activity1)
            else:
                same_assay=True
                act1 = [entry["value"] for entry in meas[pair[0]] if entry["assay_id"] in overlap]
                act2 = [entry["value"] for entry in meas[pair[1]] if entry["assay_id"] in overlap]
                activities1.append(sum(-math.log10(i/1E9) for i in act1)/len(act1))
                activities2.append(sum(-math.log10(i/1E9) for i in act2)/len(act2))
                deltas.append(activities2[-1]-activities1[-1])
            identical_assay.append(same_assay)
        transforms[tf] = dict(ligand_ids = transforms[tf], activities1 = activities1, activities2 = activities2, assay_identity = identical_assay, deltas = deltas)

# Remove pairs with identical diffs from the same assay (Database errors)
    if not keep_id_diffs:
        for tf in transforms.keys():
            if len(transforms[tf]["deltas"]) == len(set(transforms[tf]["deltas"])): continue     # Everything alright here
            act_counts = Counter(transforms[tf]["deltas"])
            new_dict = dict(ligand_ids=["" for i in act_counts], activities1=[0 for i in act_counts], activities2 = [0 for i in act_counts], assay_identity = [False for i in act_counts], deltas = act_counts.keys())
            for i in range(len(transforms[tf]["deltas"])):
                idx = new_dict["deltas"].index(transforms[tf]["deltas"][i])
                if new_dict["assay_identity"][idx] == True: continue
                new_dict["ligand_ids"][idx] = transforms[tf]["ligand_ids"][i]
                new_dict["activities1"][idx] = transforms[tf]["activities1"][i]
                new_dict["activities2"][idx] = transforms[tf]["activities2"][i]
                new_dict["assay_identity"][idx] = transforms[tf]["assay_identity"][i]
            transforms[tf] = new_dict

    return transforms


def join_transforms(transforms_Ki,transforms_IC50):
    """
Join Transformation Information from Ki and IC50 datasets
    """
    joint_transforms = transforms_Ki
    for transf,pairs in transforms_IC50.iteritems():
        if not transf in transforms_Ki.keys():
            joint_transforms[transf] = pairs
            continue
        else:
            joint_transforms[transf]["ligand_ids"] = joint_transforms[transf]["ligand_ids"] + pairs["ligand_ids"]
            joint_transforms[transf]["activities1"] = joint_transforms[transf]["activities1"] + pairs["activities1"]
            joint_transforms[transf]["activities2"] = joint_transforms[transf]["activities2"] + pairs["activities2"]
            joint_transforms[transf]["assay_identity"] = joint_transforms[transf]["assay_identity"] + pairs["assay_identity"]
            joint_transforms[transf]["deltas"] = joint_transforms[transf]["deltas"] + pairs["deltas"]

# Remove pairs with identical diffs
    for tf in joint_transforms.keys():
        if len(joint_transforms[tf]["deltas"]) == len(set(joint_transforms[tf]["deltas"])): continue     # Everything alright here
        act_counts = Counter(joint_transforms[tf]["deltas"])
        new_dict = dict(ligand_ids=["" for i in act_counts], activities1=[0 for i in act_counts], activities2 = [0 for i in act_counts], assay_identity = [False for i in act_counts], deltas = act_counts.keys())
        for i in range(len(joint_transforms[tf]["deltas"])):
            idx = new_dict["deltas"].index(joint_transforms[tf]["deltas"][i])
            if new_dict["assay_identity"][idx] == True: continue
            new_dict["ligand_ids"][idx] = joint_transforms[tf]["ligand_ids"][i]
            new_dict["activities1"][idx] = joint_transforms[tf]["activities1"][i]
            new_dict["activities2"][idx] = joint_transforms[tf]["activities2"][i]
            new_dict["assay_identity"][idx] = joint_transforms[tf]["assay_identity"][i]
        joint_transforms[tf] = new_dict

    return joint_transforms


def write_transforms_to_file(transforms,filename="dummy_transforms.txt",min_pairs=3,p_level=0.05,std_min=0.0,id_assays=True,full_info=False):
    """
Write selected transformations to file.
min_pairs  : Minimum number of pairs per transformations
p_level    : Maximum p_value
std_min    : Minimum Standard deviation of differences within pairs
id_assays  : separately output statistics for using pairs from identical assays only
    """

    print "Writing significant transformations to file"
    if min_pairs < 2:
        print "At least 2 pairs per transformation are necessary for significance tests."
        print "min_pairs set to 2"
        min_pairs = 2

    header = "Transformation\tAssay_specific\tp-value\tAverage_Activity_Difference\tSigma_Differences\tnpairs"
    if full_info: header = header+"\tLigand_IDs\tlog(Activities[nM])\tAssay_Identity"
    header = header+"\n"
    f = open(filename,"w")
    f.write(header)

    for transf,pairs in transforms.iteritems():
        if len(pairs["ligand_ids"]) < min_pairs: continue
        diffs = pairs["deltas"]
        npairs_all = len(diffs)
        p_all = stats.ttest_rel(diffs,[0.0 for i in diffs])[1]
        av_all = sum(diffs)/npairs_all
        std_all = stats.tstd(diffs)
        if npairs_all >= min_pairs and p_all <= p_level and std_all >= std_min:
            f.write(transf+"\t"+"mixed_assays"+"\t"+"{:4.2}".format(p_all)+"\t"+"{:4.3}".format(av_all)+"\t"+"{:4.2}".format(std_all)+"\t"+str(npairs_all))
            if full_info:
                for i in range(npairs_all): f.write("\t"+pairs["ligand_ids"][i][0]+":"+pairs["ligand_ids"][i][1])
                for i in range(npairs_all): f.write("\t"+"{:4.3}".format(pairs["activities1"][i])+":"+"{:4.3}".format(pairs["activities2"][i]))
                for i in range(npairs_all): f.write("\t"+str(pairs["assay_identity"][i]))
            f.write("\n")
        if id_assays == False: continue
        diffs_id = list(set([pairs["deltas"][i] for i in range(npairs_all) if pairs["assay_identity"][i]]))
        npairs_id = len(diffs_id)
        if npairs_id < min_pairs: continue
        p_id = stats.ttest_rel(diffs_id,[0.0 for i in diffs_id])[1]
        av_id = sum(diffs_id)/npairs_id
        std_id = stats.tstd(diffs_id)
        if npairs_id >= min_pairs and p_id <= p_level and std_id >= std_min:
            f.write(transf+"\t"+"ident_assays"+"\t"+"{:4.2}".format(p_id)+"\t"+"{:4.3}".format(av_id)+"\t"+"{:4.2}".format(std_id)+"\t"+str(npairs_id))
            if full_info:
                for i in range(npairs_all):
                    if pairs["assay_identity"][i] == True: f.write("\t"+pairs["ligand_ids"][i][0]+":"+pairs["ligand_ids"][i][1])
                for i in range(npairs_all):
                    if pairs["assay_identity"][i] == True:f.write("\t"+"{:4.3}".format(pairs["activities1"][i])+":"+"{:4.2}".format(pairs["activities2"][i]))
                for i in range(npairs_all):
                   if pairs["assay_identity"][i] == True:f.write("\t"+str(pairs["assay_identity"][i]))
            f.write("\n")

    f.close()


def build_ligand_dictionary_from_infile(infile):
    """
    This routine assumes that the activity is given as pActivity (=-log10(Activity [M]))
    """
    f = open(infile,'r')
    data = f.readlines()
    f.close()
  
    smi_col = 1
    cpd_id = 0
    act_col = 2
    assay_id = ''

    header = data[0].split()
    for i,item in enumerate(header):
        if "smiles" in item: smi_col = i
        if "act" in item or "log" in item: act_col = i
        if "ID" in item or "name" in item: cpd_id = i
        if "assay" in item: assay_id = i

    data = data[1:]
    meas = dict()
    for line in data:
        line = line.split()
        if line[cpd_id] in meas.keys():
            meas[line[cpd_id]].append(dict(assay_id=(line[assay_id] if assay_id != "" else "unknown"),value=10**((-1)*float(line[act_col])+9),smiles=line[smi_col].replace("\\\\","\\")))
        else:
            meas[line[cpd_id]]= [dict(assay_id=(line[assay_id] if assay_id != "" else "unknown"),value=10**((-1)*float(line[act_col])+9),smiles=line[smi_col].replace("\\\\","\\"))]

    return meas


def main():
    """
Main program routine calling all subroutines
    """
    import sys, os.path

# Defaults
    target = ""
    mysql_user = "root"
    mysql_password = ""
    use_old_results = False
    activity_type = "Ki"
    get_all_targets = False
    min_activities = 2000
    max_diff = 2.5
    outfile = ""
    min_pairs = 3
    p_level = 0.05
    std_min = 0.0
    id_assays = False
    full_info = False
    max_heavy = 70
    min_conf = 4
    keep_id_diffs = False
    mmpa_home="/usr/local/share/RDKit/Contrib/mmpa/"  # Change this accordingly if you are not Christian Kramer
    smiin = ""

    Usage = """
    Usage
    python MMP_analysis.py [-target id] [-mysql_user id] [-mysql_password id] [-use_old_results]
          [-activity_type {Ki, IC50, both}] [-get_all_targets] [-min_Nactivities #] [-max_diff #]
          [-out filename] [-npairs #] [-p_level #] [-std_min #] [-id_assays] [-full_info]
          [-max_heavy #] [-min_conf #] [-smiin file] [-keep_id_diffs]
    
    Keyword           Default         Meaning
    target            none            CHEMBL Target ID
    mysql_user        root            MYSQL User ID
    mysql_password    none            MYSQL User Password
    use_old_results   False           Use raw MMP results from previous computation
    activity_type     Ki              Select activity type to be used from [Ki, IC50, both]
    get_all_targets   False           Fetch all target IDs from CHEMBL with specified activity type
    min_Nactivities   3000            Minimum number of entries for get_all_targets
    max_diff          2.5             Maximum tolerated difference in log(Activity) for multiple
                                      measurements of a given ligand/protein complex.
    npairs            3               Minimum number of pairs per transformation
    p_level           0.05            Significance level for outputting transformations found
    std_min           0.0             Minimum standard deviation of activity differences
    id_assays         False           Additionally write results for pairs from the same assay only
    full_info         False           Output ligands IDs, activities and assay identity
    max_heavy         70              Maximum number of heavy atoms per ligand
    min_conf          4               Minimum CHEMBL Confidence Score
    smiin             none            Read smiles file instead of extracting data from ChEMBL
    keep_id_diffs     False           Do not filter out pairs which have exactly the same diff (could be database errors)

    out               target_activity_transformations.txt          Output file name
    
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
        elif sys.argv[i] == "-get_all_targets":
            get_all_targets = True
        elif sys.argv[i] == "-min_Nactivities":
            min_activities = int(sys.argv[i+1])
        elif sys.argv[i] == "-max_heavy":
            max_heavy = int(sys.argv[i+1])
        elif sys.argv[i] == "-max_diff":
            max_diff = float(sys.argv[i+1])
        elif sys.argv[i] == "-out":
            outfile = sys.argv[i+1]
        elif sys.argv[i] == "-npairs":
            min_pairs = int(sys.argv[i+1])
        elif sys.argv[i] == "-p_level":
            p_level = float(sys.argv[i+1])
        elif sys.argv[i] == "-std_min":
            std_min = float(sys.argv[i+1])
        elif sys.argv[i] == "-min_conf":
            min_conf = int(sys.argv[i+1])
        elif sys.argv[i] == "-id_assays":
            id_assays = True
        elif sys.argv[i] == "-full_info":
            full_info = True
        elif sys.argv[i] == "-keep_id_diffs":
            keep_id_diffs = True
        elif sys.argv[i] == "-smiin":
            smiin = sys.argv[i+1]
        elif sys.argv[i] == "-h":
            print Usage
            exit(0)


    if len(sys.argv) < 2:
        print Usage
        exit(0)

    if target == "" and not get_all_targets:
        print "Target not recognized. Exiting"
        print Usage
        exit(0)

    if activity_type not in ["Ki","IC50","both","other"]:
        print "Allowed activity types: Ki, IC50, both, other(only for SMILES input file)"
        print "Exiting"
        exit(0)

    if outfile == "": outfile = target+"_"+activity_type+"_transformations.txt"

    if get_all_targets:
        data = get_all_chembl_activities()
        most_frequent_targets = get_most_frequent_targets(data,min_activities)
        for id in most_frequent_targets: print id
        exit(0)

    if smiin != "":
        meas = build_ligand_dictionary_from_infile(smiin)
        if use_old_results:
            import os.path
            if not os.path.exists(target+"_"+activity_type+"_mmp_raw.csv"): use_old_results = False

        if not use_old_results:
            write_smiles_id_file(meas,target,activity_type,max_heavy)
            calc_raw_MMPs(target,activity_type,mmpa_home="/usr/local/share/RDKit/Contrib/mmpa/")

        mmps = read_raw_mmps(target,activity_type)

        transforms = workup_raw_transformations(mmps,meas,keep_id_diffs)
        write_transforms_to_file(transforms,filename=outfile,min_pairs=min_pairs,p_level=p_level,std_min=std_min,id_assays=id_assays,full_info=full_info)

        exit(0)

    elif activity_type == "both":
# Ki data first
        sdata = extract_all_activities_for_target(target,activity_type="Ki",user=mysql_user,password=mysql_password)
        sdata = clean_chembl_data(sdata, min_conf=min_conf)
        meas = build_ligand_dictionary(sdata, maxdiff=max_diff)

        if not use_old_results or not os.path.exists(target+"_Ki_mmp_raw.csv"):
            write_smiles_id_file(meas,target,"Ki",max_heavy)
            calc_raw_MMPs(target,"Ki",mmpa_home="/usr/local/share/RDKit/Contrib/mmpa/")

        mmps = read_raw_mmps(target,"Ki")
        transforms_Ki = workup_raw_transformations(mmps,meas,keep_id_diffs)
# IC50 data next
        sdata = extract_all_activities_for_target(target,activity_type="IC50",user=mysql_user,password=mysql_password)
        sdata = clean_chembl_data(sdata, min_conf=min_conf)
        meas = build_ligand_dictionary(sdata, maxdiff=max_diff)

        if not use_old_results or not os.path.exists(target+"_IC50_mmp_raw.csv"):
            write_smiles_id_file(meas,target,"IC50",max_heavy)
            calc_raw_MMPs(target,"IC50",mmpa_home="/usr/local/share/RDKit/Contrib/mmpa/")

        mmps = read_raw_mmps(target,"IC50")
        transforms_IC50 = workup_raw_transformations(mmps,meas,keep_id_diffs)
# Join Information and write to output
        joint_transforms = join_transforms(transforms_Ki,transforms_IC50)
        write_transforms_to_file(joint_transforms,filename=outfile,min_pairs=min_pairs,p_level=p_level,std_min=std_min,id_assays=id_assays,full_info=full_info)

        exit(0)


    else:
        sdata = extract_all_activities_for_target(target,activity_type=activity_type,user=mysql_user,password=mysql_password)
        sdata = clean_chembl_data(sdata, min_conf=min_conf)
        meas = build_ligand_dictionary(sdata, maxdiff=max_diff)

        if use_old_results:
            import os.path
            if not os.path.exists(target+"_"+activity_type+"_mmp_raw.csv"): use_old_results = False

        if not use_old_results:
            write_smiles_id_file(meas,target,activity_type,max_heavy)
            calc_raw_MMPs(target,activity_type,mmpa_home="/usr/local/share/RDKit/Contrib/mmpa/")

        mmps = read_raw_mmps(target,activity_type)

        transforms = workup_raw_transformations(mmps,meas,keep_id_diffs)
        write_transforms_to_file(transforms,filename=outfile,min_pairs=min_pairs,p_level=p_level,std_min=std_min,id_assays=id_assays,full_info=full_info)

        exit(0)


if __name__ == "__main__":
  main()



