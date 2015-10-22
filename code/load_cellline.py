import pandas as pd
import numpy as np
import pickle
import os
import pdb
import re

def load_gdsc(cell_lines_dir):
    gdsc = pd.read_table(cell_lines_dir + 'cancerrxgene/gdsc_drug_sensitivity_fitted_data_w5.csv', sep=',',na_values='null')
    manov = pd.read_table(cell_lines_dir + 'cancerrxgene/gdsc_manova_output_w5.csv', 
                        sep=',',usecols=['drug','drug_id']).drop_duplicates()
    gdsc = pd.merge(gdsc, manov,on='drug_id',how='left')
    gdsc['drug'] = gdsc['drug'].str.upper()
    return gdsc

def load_ctrp(cell_lines_dir):
    ctrp = pd.read_table(cell_lines_dir + 'CTRP/tableS2G_sensitivity_AUC.txt',sep='\t')
    ctrp['drug'] = ctrp['compound_name'].str.upper()
    return ctrp

## load data frame of synonyms 
def synonym_load(synonym_pkl):
    resp_to_syn = pickle.load(open(synonym_pkl))
    #import pdb
    #pdb.set_trace()
    dict_of_syn = dict()
    for to_syn in resp_to_syn:
        ## first, remove extra hits try to keep it to a main hit
        for d in resp_to_syn[to_syn]:
            if len(resp_to_syn[to_syn][d]) > 1:
                if d in resp_to_syn[to_syn][d]:
                   resp_to_syn[to_syn][d] = set([d])
                else:
                    rerep_d = re.sub(r"[\W_]", "", d)
                    rerep_s = set([re.sub(r"[\W_]", "", s)
                                for s in resp_to_syn[to_syn][d]])  
                    if rerep_d in rerep_s:
                        resp_to_syn[to_syn][d] = set([rerep_d])
        hits = pd.DataFrame(dict([(d, list(resp_to_syn[to_syn][d])) 
                                    for d in resp_to_syn[to_syn]
                                    if len(resp_to_syn[to_syn][d]) == 1]),
                            index=['syn']).transpose()
        #pdb.set_trace()
        for d in resp_to_syn[to_syn][d]:
            if len(resp_to_syn[to_syn][d]) > 1:
                hits = pd.concat((hits,
                                 pd.DataFrame(dict([(d, l) for l in list(resp_to_syn[to_syn][d])]), index=['syn'].tranpose())), axis=0)
        dict_of_syn[to_syn] = hits
    return dict_of_syn

def synonym_prepare(dict_of_df, target_list, pkl_name):
    dict_of_syn_lists = dict([(d, dict_of_df[d]['drug'].unique())
                               for d in dict_of_df])
    target = set(target_list) # set(target['drug'])
    pkl = open(pkl_name + '.pkl', 'w')
    pickle.dump((dict_of_syn_lists, target), pkl)
    pkl.close()

def synonym_add(dict_of_df, dict_of_syn):
    dict_of_syn_added = dict()
    for to_add in dict_of_df:
        dict_of_syn_added[to_add] = pd.merge(dict_of_df[to_add],
                                             dict_of_syn[to_add],
                                             left_on='drug',
                                             right_index=True,how='left')
    return dict_of_syn_added
