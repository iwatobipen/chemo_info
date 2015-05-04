import pandas as pd
import numpy as np
import sys, os

#SAVE the smiles for RDKIT_MMPA.
data = pd.read_table( "ks_compound.txt", sep = "\t")
data = pd.DataFrame( data = data, columns= ["SMILES","COMPOUND_ID"] )
data.to_csv( "chembl.smi", sep = "\t", index = False, cols = [ "SMILES", "COMPOUND_ID" ] )

#read activity data and data clearning
data = pd.read_table( "ks_bioactivity.txt", sep = "\t" )
data = data[data.ACTIVITY_TYPE == "IC50"]
data = data[data.ASSAY_TYPE == "B"]
data = data[data.RELATION == "="]
data = pd.DataFrame(data = data, columns =[ "ACTIVITY_ID", "NAME", "COMPOUND_ID","ACTIVITY_TYPE","STANDARD_VALUE", "STANDARD_UNIT"])
#data.to_csv( "datadata.txt" , index=False)

#########################################
#separate activity data by each targets.
#########################################
target_set = []
for target in data.NAME:
    if target in target_set:
        pass
    else:
        target_set.append(target)

print len(target_set)

##########################################
# make target vs activity data dictionaly.
##########################################
splitdata = {}

for target in target_set:
    splitdata[target] = data[data.NAME == target]

#########################################
#make mmpa_set
#########################################
os.system( "python rfrag.py < chembl_nohead.smi > chembl_frag.txt" )
os.system( "python index.py < chembl_frag.txt > chembl_mmpa.txt" )
#now get mmps "chembl_mmpa.txt"
mmpa_data = pd.read_table("chembl_mmpa.txt", sep = "," ,header=None, names = ["SMILES_L","SMILES_R","ID_L","ID_R","TRANSFORM", "CONTEXT" ])

# make mmp_data frame that has each target activity data.
mmpa_act_dict = {}

for k, v in splitdata.items():
    merge_data = mmpa_data.merge( v, left_on = "ID_L", right_on ="COMPOUND_ID",how = "inner")
    merge_data = merge_data.rename( columns={ "ACTIVITY_ID":"ACTIVITY_ID_L",
                                   "NAME":"NAME_L",
                                   "COMPOUND_ID":"COMPOUND_ID_L",
                                   "ACTIVITY_TYPE":"ACTIVITY_TYPE_L",
                                   "STANDARD_VALUE":"STANDARD_VALUE_L",
                                   "STANDARD_UNIT":"STANDARD_UNIT_L"})

    merge_data = merge_data.merge( v, left_on = "ID_R", right_on="COMPOUND_ID",how = "inner")
    merge_data = merge_data.rename( columns={ "ACTIVITY_ID":"ACTIVITY_ID_R",
                                   "NAME":"NAME_R",
                                   "COMPOUND_ID":"COMPOUND_ID_R",
                                   "ACTIVITY_TYPE":"ACTIVITY_TYPE_R",
                                   "STANDARD_VALUE":"STANDARD_VALUE_R",
                                   "STANDARD_UNIT":"STANDARD_UNIT_R"})
    mmpa_act_dict[ k ] = merge_data


df_list = mmpa_act_dict.values()
final_data = pd.concat(df_list)
final_data.to_csv( "result.txt", index=False )


