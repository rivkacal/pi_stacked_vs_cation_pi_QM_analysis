'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Details:
For a given residue_i_name and residue_j_name look for all pairs within some distance, apply rotation and translation
to find the following angles where values taken from ****:
Fully stacked defenition: P < 30, Ttheta1 > 60 Ttheta2 > 60
Staggered stacking:       P < 30, 30 < Ttheta1,Ttheta2 < 60
Parallel in-plane:        P < 30, Ttheta1,Ttheta2 < 30
Tilted:              30 < P < 70,
Edge-ring-face            P > 70, Ttheta1 > 60, Ttetha2 > 30
Cogwheel                  P > 70, 0 < Ttheta1,Ttheta2 < 60

*****for parameters please see PHE-PHE geometry search paper:
https://reader.elsevier.com/reader/sd/pii/0014579385809820?token=E75C9BBEDB5EDF547E41B1FD6290C4C1FA8CB25F59B0E8B5DAD7043060819E1168B7A1A691B41CE93F1ECF7710455EFC&originRegion=eu-west-1&originCreation=20221117162337

Use:
python specify_res_i_res_j_geometry_phi.py <res_i_name> <res_j_name> <csv_out_name.csv>

Example:
python specify_res_i_res_j_geometry_phi.py HIS PHE HIS_PHE_params_phi.csv
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
res_i_name = sys.argv[1]
res_j_name = sys.argv[2]
csv_out_name = sys.argv[3]
# name for output csv file saving pdb files name that contain a matching pair(s) and the pairs ids (serial num)
Point = namedtuple("Point", ["x", "y"])

global dict_df  # dataframe that updates with matching pairs ids above
global threshold  # threshold for closest atom distance (perliminary) scan for pairs

shortest_thresh = 4.6  # as in PHE-PHE paper
centroid_dist_thresh = 10

############################################ END OF SECTION #######################################################################################################################

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''


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

def specify_X_Y_given_geomtery(res_i_name, res_j_name, chains_dict, shortest_thresh, centroid_dist_thresh, title,
                               list_pairs, struct, pdbs_path):

    reference_res_i, reference_coord_i = pdc.assign_reference_res(res_i_name, pdbs_path)
    reference_res_j, reference_coord_j = pdc.assign_reference_res(res_j_name, pdbs_path)
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

                '''~~~~~~CHECK res_j_name~~~~~~~~~~~~'''

                # choose only interesting res j later
                for residue_j in res_j_arr:
                    residue_j_idx = residue_j.get_id()[1]
                    if abs(residue_j_idx - his_i_idx) > 0:

                        if res_j_name == 'HIS':
                            # if the second residue is also HIS then can run only over j < i and not twice ...
                            if residue_j_idx in histag_indices or residue_j_idx <= his_i_idx or residue_j_idx in list_metal_bound_idx :
                                continue
                                # print(dict_df)
                        list_pairs = pdc.assign_X_Y_phi_geometry_params_Hb(his_i, his_i_idx, residue_j, residue_j_idx,
                                                             shortest_thresh, centroid_dist_thresh, title, chain_id,
                                                             list_pairs,reference_res_i, reference_res_j)
        else:
            '''~~~~~~CASE B: res_i_name != 'HIS' ~~~~~~~~~~~~'''
            # run over all histidines i and not all residues due to unnecessary computational cost
            for residue_i in res_i_arr:
                # for every histidine i get the residue index for hydrophobic calculations and to check if it is in histag group to ignore
                residue_i_idx = residue_i.get_id()[1]

                '''~~~~~~CHECK res_j_name~~~~~~~~~~~~'''
                # only for histidines look for histags to remove...

                for residue_j in res_j_arr:
                    residue_j_idx = residue_j.get_id()[1]
                    if res_i_name == res_j_name:  # both same residues but not HIS, so donot count twice:
                        if residue_j_idx - residue_i_idx > 0:
                            list_pairs = pdc.assign_X_Y_phi_geometry_params_Hb(residue_i, residue_i_idx, residue_j, residue_j_idx,
                                                                 shortest_thresh, centroid_dist_thresh, title, chain_id,
                                                                 list_pairs,reference_res_i, reference_res_j)

                    elif abs(residue_j_idx - residue_i_idx) > 0:

                        if res_j_name == 'HIS':
                            histag_indices = lst.check_histag(res_j_arr)
                            if residue_j_idx in histag_indices:
                                continue
                            if find_nearby_cations(residue_j):
                                list_metal_bound.append(residue_j)
                                list_metal_bound_idx.append(residue_j_idx)
                                continue

                        list_pairs = pdc.assign_X_Y_phi_geometry_params_Hb(residue_i, residue_i_idx, residue_j, residue_j_idx,
                                                             shortest_thresh, centroid_dist_thresh, title, chain_id,
                                                             list_pairs,reference_res_i, reference_res_j)

    return list_pairs

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ANALYSIS END~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
###################################################################################################################################################################################################

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

# some status variables
file_counter = 0
list_pairs = []

# scanning all files in full path provided for matching pairs:
for file in os.listdir(full_path):
    if file.endswith('.pdb.gz'):
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
        list_pairs = specify_X_Y_given_geomtery(res_i_name,res_j_name,chains_dict, shortest_thresh, centroid_dist_thresh, title, list_pairs, struct, full_path)

        os.remove(pdb_file)
        file_counter += 1
        print(file_counter)
        cprint(f"{pdb_file} finished", "blue")

dict_df = pd.DataFrame(list_pairs)
dict_df.to_csv(csv_out_name)
