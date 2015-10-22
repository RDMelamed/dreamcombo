import time
import pickle
import pandas as pd
import numpy as np
from itertools import chain
import pdb

def match_drug_resp_lincs(pkl_file,
                          inst_info_file):
    (df_dict, # dict of drug response source -> df of [drug,cell,auc]
    cp_df, # DF of perturbations
    unified_2_lincs_pertname, # mapping of unified drug -> lincs drug
    unified_2_lincs_cell # unified cell -> lincs cell
    ) = pickle.load(open(pkl_file))

    ## get common set of drugs
    drugs_shared = set(unified_2_lincs_pertname.keys()) & \
      set(chain.from_iterable([df_dict[df]['drug'] for df in df_dict]))
    cells_shared = set(unified_2_lincs_cell.keys()) & \
      set(chain.from_iterable([df_dict[df]['cell_line_name'] for df in df_dict]))
    lincs_cell_2_shared_unified = dict([(unified_2_lincs_cell[c], c)
                                        for c in cells_shared])
    rep = open('rep.drug.txt','w')
                        
    t0 = time.time()
    # get list of perturbations with shared drug & cells
    for df in df_dict:
        df_dict[df] = df_dict[df].rename(columns={'cell_line_name':'unified_cell',
                                            'drug':'unified_drug',
                                            'auc':df + 'auc'})

    ## get pert name -> unified drug, unified cellS
    pert_to_drug_cells = dict()
    for drug in drugs_shared:
        drug_df = pd.DataFrame() #np.zeros())
        pert_iname_mat = cp_df.loc[cp_df['pert_iname']==unified_2_lincs_pertname[drug],:]
        for id_row in range(pert_iname_mat.shape[0]):
            got_cell_matches = pert_iname_mat.iloc[id_row,:]['cell_id']
            pert = pert_iname_mat.iloc[id_row,:]['pert_id']
            for shared_cell in [c for c in cells_shared
                                if unified_2_lincs_cell[c] in got_cell_matches]:
                if not pert in pert_to_drug_cells:
                    pert_to_drug_cells[pert] = [[drug], set()]
                if not drug in pert_to_drug_cells[pert][0]:
                    rep.write('drug-msimatch new_drug=' + drug + ' old_drug=' + ','.join(pert_to_drug_cells[pert][0]) + ' perturb=' + pert + '\n')
                    pert_to_drug_cells[pert][0].append(drug)
                pert_to_drug_cells[pert][1] |= set([shared_cell])


    # first create df of <pert_id> <unified_drug>
    # then for each one, add in a row for each <unified_cell>
    # as well as <cell_id> ( this is the lincs cell name)
    pert2 = pd.DataFrame(dict([(p, list(pert_to_drug_cells[p])[0])
                               for p in pert_to_drug_cells])).transpose()
    pert2['pert'] = pert2.index
    pertlist = []
    for pert in pert2.index:
        cells = list(pert_to_drug_cells[pert][1])
        pert_drug = pd.DataFrame(np.tile(pert2.loc[pert,:], (len(cells),1)),
                                 columns=['unified_drug','pert_id'])
        pert_drug['unified_cell'] = cells
        pert_drug['cell_id'] = [unified_2_lincs_cell[c] for c in cells]
        
        pertlist.append(pert_drug)
    pert_to_drug_cells = pd.concat(pertlist, axis=0)
    
    ## add in DMSO -- still want it though is not a pert
    dmso_mat = pert_to_drug_cells.loc[:,['unified_cell','cell_id']].drop_duplicates()
    dmso_mat['unified_drug'] = 'DMSO'
    dmso_mat['pert_id'] = 'DMSO'
    pert_to_drug_cells = pd.concat((pert_to_drug_cells, dmso_mat), axis=0)

    header = open(inst_info_file).readline().strip().split('\t')
    added_col = ['unified_drug','unified_cell'] + [df + 'auc' for df in df_dict]
    ## prepare output files
    for drug in drugs_shared | set(['DMSO']):
        with open ('drug_inst/' + drug,'w') as f:
            #drug_res_inst[drug] = open
            f.write('\t'.join(header + added_col) + '\n')

    for cell in cells_shared:
        with open ('cell_inst/' + cell,'w') as f:
            #drug_res_inst[drug] = open
            f.write('\t'.join(header + added_col) + '\n')
            
    for inst_chunk in pd.read_table(inst_info_file,sep='\t',
                                   skiprows=1,index_col=None,
                                   chunksize=10000, header=None,
                                   names=header):
        ## merge in by cell & pert
        inst_chunk = pd.merge(inst_chunk,
                              pert_to_drug_cells,
                              on = ['pert_id','cell_id'],how='inner')

        ## then merge data from CTRP/GDSC
        for df in df_dict:
            inst_chunk = pd.merge(inst_chunk,
                                   df_dict[df],
                                   on = ['unified_drug','unified_cell'],
                                   how='left')

        ## now keep only those with relevant info
        keep = inst_chunk['pert_id']=='DMSO'
        for df in df_dict:
            keep = keep | ~pd.isnull(inst_chunk[df + 'auc'])
        #inst_chunk.loc[keep,'pert_id'].groupby(['pert_id']).count()
        inst_chunk = inst_chunk.loc[keep,:]

        ## now write out info from this chunk
        for drug in inst_chunk['unified_drug'].unique():
            with open('drug_inst/' + drug,'a') as f:
                inst_chunk.loc[inst_chunk['unified_drug']==drug,
                               header + added_col].to_csv(f, header=False, index=False,sep='\t',na_rep='NA')
        for cell in inst_chunk['unified_cell'].unique():
            with open('cell_inst/' + cell,'a') as f:
                inst_chunk.loc[inst_chunk['unified_cell']==cell,
                               header + added_col].to_csv(f, header=False, index=False, sep='\t',na_rep='NA')

    t1 = time.time()
    rep.write('finished in {:2.2f} minutes'.format((t1 - t0)/60))
    rep.close()
    return 

if __name__ == "__main__":
    import sys
    match_drug_resp_lincs(sys.argv[0], sys.argv[1])
