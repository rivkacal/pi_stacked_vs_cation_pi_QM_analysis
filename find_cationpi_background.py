'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Details:
non-HIS cation-pi background (xyz files as in find_HIS_xyz_backgorund)

Use:
python find_cationpi_background.py TRP 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
import sys
import os
from collections import namedtuple
import numpy as np
import gzip
import shutil
from Bio.PDB import *
from termcolor import cprint
import pandas as pd
# local modules:
import PlanesDistCalc as pdc
import LoadStruct as lst


############################################## GLOBALS ############################################################################################################################
full_path = '/home_d/rivka/Bioinformatics_HIS/python_codes/'  # for searching pdb files in
#full_path = '/home_d/rivka/Bioinformatics_HIS/python_codes/test_analysis/scripts/'

res_i_name = sys.argv[1]

# name for output csv file saving pdb files name that contain a matching pair(s) and the pairs ids (serial num)
Point = namedtuple("Point", ["x", "y"])

global dict_df  # dataframe that updates with matching pairs ids above

shortest_thresh = 6  # as in PHE-PHE paper

############################################ END OF SECTION #######################################################################################################################

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''


#metals screening for HIS only (ligand)
def find_nearby_cations(residue, ns) -> bool:
    if residue.has_id('CG') and residue['CG'].coord.shape == (3,):
        close_atoms = ns.search(residue['CG'].coord, 5)
        for atom in close_atoms:
            if atom.get_parent() == None:
                continue
            if atom.name == 'FE' or atom.name == 'ZN' or atom.name == 'CU' or (
                    atom.name == 'CA' and atom.get_parent().resname == 'CA') or atom.name == 'MG' or atom.name == 'NA' or atom.name == 'NI' or atom.name == 'AG' or atom.name == 'CD' or atom.name == 'CO' or atom.name == 'MN':
                # print(close_atoms)
                return True
    return False

def assign_HIS_cationpi_pairs(chains_dict, shortest_thresh, title,
                               list_pairs, struct, res_i_name, res_j_name, atom_name, is_CH_arr, is_Hb_arr, is_cationpi_arr):
    pdbs_path = '/home_d/rivka/Bioinformatics_HIS/python_codes/'

    chains_res_dict = {chain_id: dict(i=[], j=[]) for chain_id, chaing_obj in chains_dict.items()}
    # first thing create an array of residue X to run over (reduce computational cost not to run over all residues...)
    for chain_id, chain_obj in chains_dict.items():
        for residue in chain_obj.get_list():
            tmp_arr = []
            if residue.resname == res_i_name:
                tmp_arr = chains_res_dict[chain_id]['i']
                tmp_arr.append(residue)
                chains_res_dict[chain_id]['i'] = tmp_arr

            if residue.resname == res_j_name:
                tmp_arr = chains_res_dict[chain_id]['j']
                tmp_arr.append(residue)
                chains_res_dict[chain_id]['j'] = tmp_arr
    #once in a given structure find the atom lists
    atoms = Selection.unfold_entities(struct, 'A')
    ns = NeighborSearch(atoms)

    for chain_id, chain_res_dict in chains_res_dict.items():
        res_i_arr = chains_res_dict[chain_id]['i']
        res_j_arr = chains_res_dict[chain_id]['j']

        '''~~~~~~CASE A: res_i_name == 'HIS' ~~~~~~~~~~~~'''
        #we'll need to screen histidines which are in the vicinity of metals (binding imposee conformation)

        list_metal_bound = []
        list_metal_bound_idx = []
        if res_i_name == 'HIS':
            # only for histidines look for histags to remove...
            histag_indices = lst.check_histag(res_i_arr)
            # run over all histidines i and not all residues due to unnecessary computational cost
            for his_i in res_i_arr:
                # for every histidine i get the residue index for hydrophobic calculations and to check if it is in histag group to ignore
                his_i_idx = his_i.get_id()[1]
                if his_i_idx in histag_indices:
                    # ignore this HIS
                    continue

                if find_nearby_cations(his_i, ns):
                    list_metal_bound.append(his_i)
                    list_metal_bound_idx.append(his_i_idx)
                    continue

                '''~~~~~~CHECK res_j~~~~~~~~~~~~'''

                # choose only interesting res j later
                for residue_j in res_j_arr:
                    residue_j_idx = residue_j.get_id()[1]
                    if abs(residue_j_idx - his_i_idx) > 0:
                        if residue_j.resname == 'HIS':
                            # if the second residue is also HIS then can run only over j < i and not twice ...
                            if residue_j_idx in histag_indices or residue_j_idx <= his_i_idx or residue_j_idx in list_metal_bound_idx :
                                continue
                        list_pairs,is_CH_arr, is_Hb_arr, is_cationpi_arr = pdc.extract_one_coord(his_i, residue_j,pdbs_path, atom_name,list_pairs,is_CH_arr, is_Hb_arr, is_cationpi_arr)

        elif res_i_name == 'PHE' or res_i_name == 'TRP' or res_i_name == 'TYR':
            for res_i in res_i_arr:
                res_i_idx = res_i.get_id()[1]
                for residue_j in res_j_arr:
                    residue_j_idx = residue_j.get_id()[1]
                    if abs(residue_j_idx - res_i_idx) > 0:

                        list_pairs,is_CH_arr, is_Hb_arr, is_cationpi_arr = pdc.extract_one_coord(res_i, residue_j,pdbs_path, atom_name,list_pairs,is_CH_arr, is_Hb_arr, is_cationpi_arr)

    return list_pairs

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ANALYSIS END~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
###################################################################################################################################################################################################

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

# some status variables
file_counter = 0
lys_NZ_coords = []
lys_CE_coords = []
arg_CZ_coords = []
arg_NE_coords = []
arg_NH1_coords = []
arg_NH2_coords = []
arg_CD_coords = []
is_CH_arr = []
is_Hb_arr = [] 
is_cationpi_arr = []

# scanning all files in full path provided for matching pairs:
for file in os.listdir(full_path):
    if not file.endswith('.pdb.gz'):
        continue
    filename = full_path + file

    # opening file to store coordinates and atoms type
    pdb_file = filename[:-3]

    with gzip.open(filename, 'rb') as f_in:
        with open(pdb_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    struct = lst.get_pdb_structure(pdb_file)
    chains_dict = lst.create_chain_dict(struct)
    if len(chains_dict) < 1:
        cprint(f"{pdb_file} does not contain any relevant chains, continue reading next pdb", "red")
        continue
    title = pdb_file[-8:]

    cprint(f"{pdb_file} started", "blue")
    lys_NZ_coords = assign_HIS_cationpi_pairs(chains_dict, shortest_thresh, title,
                            lys_NZ_coords, struct, res_i_name, 'LYS', 'NZ',is_CH_arr, is_Hb_arr, is_cationpi_arr)
    lys_CE_coords = assign_HIS_cationpi_pairs(chains_dict, shortest_thresh, title,
                            lys_CE_coords, struct, res_i_name, 'LYS', 'CE',is_CH_arr, is_Hb_arr, is_cationpi_arr)
    arg_CZ_coords = assign_HIS_cationpi_pairs(chains_dict, shortest_thresh, title,
                            arg_CZ_coords, struct, res_i_name, 'ARG', 'CZ',is_CH_arr, is_Hb_arr, is_cationpi_arr)
    
    arg_NE_coords = assign_HIS_cationpi_pairs(chains_dict, shortest_thresh, title,
                            arg_NE_coords, struct, res_i_name, 'ARG', 'NE',is_CH_arr, is_Hb_arr, is_cationpi_arr)
    arg_NH1_coords = assign_HIS_cationpi_pairs(chains_dict, shortest_thresh, title,
                            arg_NH1_coords, struct, res_i_name, 'ARG', 'NH1',is_CH_arr, is_Hb_arr, is_cationpi_arr)
 
    arg_CD_coords = assign_HIS_cationpi_pairs(chains_dict, shortest_thresh, title,
                            arg_CD_coords, struct, res_i_name, 'ARG', 'CD',is_CH_arr, is_Hb_arr, is_cationpi_arr)
    os.remove(pdb_file)
    file_counter += 1
    print(file_counter)
    


# save these results
last_idx = 0
#all_coords_arr = [lys_NZ_coords, lys_CE_coords, arg_CZ_coords, arg_NE_coords, arg_NH1_coords, arg_NH2_coords]
all_coords_arr = [arg_CD_coords]
#atoms_arr = ['NZ', 'CE', 'CZ', 'NE', 'NH1', 'NH2']
atoms_arr = ['CD']
#cation_name = ['LYS','LYS','ARG','ARG','ARG','ARG','ARG']
cation_name = 'ARG'
for idx, coords_arr in enumerate(all_coords_arr):
    cuur_coord_len = len(coords_arr)
    data = {'X': [x[0] for x in coords_arr],
            'Y':  [x[1] for x in coords_arr],
            'Z': [x[2] for x in coords_arr],
            'is_CH': is_CH_arr[last_idx:last_idx+cuur_coord_len],
            'is_Hb': is_Hb_arr[last_idx:last_idx+cuur_coord_len],
            'is_cationpi': is_cationpi_arr[last_idx:last_idx+cuur_coord_len]}
    #print(cuur_coord_len, len(is_CH_arr[last_idx:last_idx+cuur_coord_len]))
    df = pd.DataFrame(data)
    df.to_csv(f'{res_i_name}_{cation_name[idx]}_{atoms_arr[idx]}_coords.csv')
    last_idx = last_idx + cuur_coord_len
