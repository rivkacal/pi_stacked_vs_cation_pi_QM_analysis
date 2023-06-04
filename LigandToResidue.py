'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
this scripts gets a pdb ligand input with non distinguished C,H,N,O atoms etc.. and assigns the correct residue
and atoms according to the pdb atom naming convention!

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

import sys
import numpy as np
from Bio.PDB import *
from collections import Counter

# input parameters
input_filename = sys.argv[1]  # "min_test_2HSE_solv.pdb"  #sys.argv[1] #'current_minima_2HSD.pdb'
dir_path = sys.argv[2]
# '/home_d/rivka/Bioinformatics_HIS/python_codes/test_analysis/pdb_files/QM/final_trj_HSE_HSE/'

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SECTION A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
"""The following functions help take a current pdb ligand file and numerically assign each atom, in order to distinguish
 them for the following use of Bio.PDB modules! These are used to create a mid outputfile: 
 'sorted_<inputfilename>.pdb' that can be later deleted """


# for later modifications create a subfunction to return numbered strings of input letters according to their idx

def get_ligands_num(full_path: str) -> int:
    ligands_curr = 0
    # full_path = dir_path + input_filename
    with open(full_path, 'r') as rfile:
        for line in rfile:
            if line.startswith('HETATM'):
                splited_line = line.split()
                if int(splited_line[4]) > ligands_curr:
                    ligands_curr = int(splited_line[4])
    return ligands_curr


# ############################################### STRUCTURE MODIFICATION ###############################################
def check_atom_in_group(atom: int, groups: dict):
    """ run over all groups and check if the atom appears in any. Note that the groups atoms are dynamically
    changing within this function.

    Parameters:
    atom (int) - the index of atom, starting at 1 (at pdb file serial number column)
    groups (dict) - a dict of all groups of atoms (as keys for each group idx)

    Returns:
        the group index to which this atom belongs (either integer or None)"""

    for group_name, group in groups.items():
        if atom in group:  # found which group this atom belongs to
            return group_name[-1]

    return None
def separate_xyz_ligands(ligands_num: int, connect_arr_list: list) -> np.ndarray:
    """ This function is called only if the number of requested ligands does not fit the extracted ones in sortedX.pdb.
    The separation is according to bonds-length (inter-molecular < 1.9A).

    Parameters:
    ligands_num (int) - the number of expected ligands
    connect_arr_list (list) - a list of arrays of connectivity information (bonded atoms)

    Returns:
        ligands_arr (np.ndarray) - of the ligands idx according to each atom index location in the pdb input"""

    if ligands_num < 2:
        raise ValueError('required at least 2 ligands expected for output')
    n_connect_arr = len(connect_arr_list)  # extract the size
    # create X empty arrays to store for each group. In this script only two ligands are considered
    groups = []

    groups = {}
    for i in range(1, ligands_num + 1):
        groups[f"group_{i}"] = []

    # for i in range(1, ligands_num+1):
    #     exec(f"group_{i} = []")
    #     groups.append(eval(f"group_{i}"))

    group_x = []  # store atom values yet not defined (does not belong to an existing group)
    # connections_num = len(connect_arr_list)  # how much different arrays we have
    # first step: move the first array to group
    temp_connect_arr_list = connect_arr_list.copy()
    for atom in temp_connect_arr_list[0]:
        groups['group_1'].append(atom)
        # exec(f"group_{1}.append({atom})")
    temp_connect_arr_list[0] = []  # remove this line now

    # for any array line perform the following
    curr_atom_identified = False
    for array_idx, array in enumerate(temp_connect_arr_list):
        if not array == []:
            for atom in array:
                curr_group_idx = check_atom_in_group(atom, groups)
                if curr_group_idx is not None:   # which means atom and all these array belongs to group i
                    curr_atom_identified = True
                    break  # no need to check other atoms in this array so assign all the atoms now:

            if curr_atom_identified:
                for atom in array:
                    if not atom in groups[f"group_{curr_group_idx}"]:
                    # exec(f"group_{curr_group_idx}.append({atom})") # append to correct group
                        groups[f"group_{curr_group_idx}"].append(atom)
                temp_connect_arr_list[array_idx] = []  # remove from the running array
                curr_atom_identified = False
            else:
                for atom in array:
                    if not atom in group_x:
                        group_x.append(atom)

    # in this function yet to implement more than 2 ligand outputs. So far all group 1 will be in group_1 array,
    # and the second group terms will be on group_x
    i = 1
    group_1 = groups[f"group_{i}"]
    group_2 = group_x

    return group_1, group_2

# ############################################ STRUCTURE READING AND NAMING ############################################
def assign_letter_index(letter: str, current_counter: str) -> str:
    return letter + current_counter

# for a given ligand id (first='1', second ='2' etc ...) return an array of replaced lines with numbered atoms:
def number_ligand_atoms(full_path: str, current_ligand_idx: str) -> list:
    ligands_num = get_ligands_num(full_path)
    single_ligand_identified = ligands_num == 1
    C_counter = 0
    N_counter = 0
    H_counter = 0
    O_counter = 0
    S_counter = 0
    lines_arr = []
    connectivity_arr = []  # store all the lines starting with connectivity
    connect_arr_arr = []  # array of arrays of connectivity
    with open(full_path, 'r') as rfile:
        for line in rfile:
            if line.startswith('HETATM'):
                splited_line = line.split()
                if splited_line[4] == current_ligand_idx:
                    if splited_line[2] == 'C':
                        C_counter += 1
                        new_line = line.replace(' C', assign_letter_index('C', str(C_counter)), 1)
                        if C_counter >= 10:
                            new_line = new_line[0:15] + new_line[16:]
                    elif splited_line[2] == 'N':
                        N_counter += 1
                        new_line = line.replace(' N', assign_letter_index('N', str(N_counter)), 1)
                        if N_counter >= 10:
                            new_line = new_line[0:15] + new_line[16:]
                    elif splited_line[2] == 'H':
                        H_counter += 1
                        new_line = line.replace(' H', assign_letter_index('H', str(H_counter)), 1)
                        if H_counter >= 10:
                            new_line = new_line[0:15] + new_line[16:]
                    elif splited_line[2] == 'O':
                        O_counter += 1
                        new_line = line.replace(' O', assign_letter_index('O', str(O_counter)), 1)
                        if O_counter >= 10:
                            new_line = new_line[0:15] + new_line[16:]
                    elif splited_line[2] == 'S':
                        new_line = line.replace(' S', assign_letter_index('S', str(S_counter)), 1)
                        S_counter +=1
                    else:
                        raise ValueError('The file contains atom that is not C,N,O,S or H')
                    lines_arr.append(new_line)

            elif line.startswith('CONECT') and single_ligand_identified:
                connectivity_arr.append(line[9:-1])

    if single_ligand_identified:
        # reorganize the data to separate arrays
        for array in connectivity_arr:
            split_arr = array.split(' ')
            int_arr = []
            for element in split_arr:
                if not element == '' and not element == ' ':
                    int_arr.append(element)
            connect_arr_arr.append(int_arr)

        required_ligands = 2  #studying here only pairwise interactions
        ligand1_ids, ligand2_ids = separate_xyz_ligands(required_ligands, connect_arr_arr)
        # print(ligand1_ids)
        # print('***')
        # print(ligand2_ids)
        final_lines_arr = []
        #again for 2 ligands only replace those who do not belong
        for line_idx, line in enumerate(lines_arr):
            atom_idx = line_idx + 1
            if str(atom_idx) in ligand1_ids:
                final_lines_arr.append(line)
                continue
            else:
                new_line = line[:20] + '2' + line[21:25] + '2' + line[26:]
                final_lines_arr.append(new_line)


    rfile.close()
    if single_ligand_identified:
        return final_lines_arr
    return lines_arr


# do this for all ligands and create the sorted output file!
def create_output(ligands_num: int, input_filename: str, dir_path):
    full_path = dir_path + input_filename
    final_arr = []
    for i in range(1, ligands_num + 1):
        temp_arr = number_ligand_atoms(full_path, str(i))
        final_arr.append(temp_arr)

    output_path = dir_path + 'sorted_' + input_filename
    with open(output_path, 'w') as ofile:
        for ligand in final_arr:
            for line in ligand:
                ofile.write(line)
    ofile.close()


'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END of SECTION A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SECTION B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
'''The following functions loads the temporarily 'sorted' outputfile using Bio.PDB modules!  '''


########################################################## STURCTURE LOADING ###########################################################
def get_ligand_size(full_path: str, current_ligand_idx: int) -> int:
    counter = 0
    with open(full_path, 'r') as rfile:
        for line in rfile:
            if line.startswith('HETATM'):
                splited_line = line.split()
                if splited_line[4] == str(current_ligand_idx):
                    counter += 1
    return counter


def get_pdb_structure(pdb_file_name: str) -> Structure.Structure:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("current_struct", pdb_file_name)
    return structure

def create_chain_dict(structure: Structure.Structure) -> dict:
    ''' creates a dictionary of protein chains from a given biopython pdb structure. The keys are
        the chain ids, the values are the Chain objects'''
    model = structure[0]
    chain_dict = {}
    for chain in model.get_list():
        chain_dict[chain.get_full_id()[2]] = chain

    return chain_dict


################################################### BONDS MATRIX CALCULATONS #############################################################

def two_pt_dist(point1, point2):
    ''' This function returns the Euclidean distance between two input points '''
    # calculating Euclidean distance
    # using linalg.norm()
    dist = np.linalg.norm(point1 - point2)
    return dist


def create_inner_dist_matrix(residue, size):
    matrix = np.eye(size)
    for idx1, atom1 in enumerate(residue):
        for idx2, atom2 in enumerate(residue):
            dist = two_pt_dist(atom1, atom2)
            matrix[idx1, idx2] = dist
    return matrix


#################################################### RESIDUE & BONDS IDENTIFICATION ######################################################

def find_bond_type(inner_dist_matrix, sorted_idx_mat, row_idx, atoms_dict):
    sigma = 0.0655  # current 'error'
    sigma_CdO = 0.9
    sigma_NH = 0.72

    # for aromatic add lowercase ar. lengths taken from: https://hydra.pharm.monash.edu/modules/mod2/bondlen.html
    # https://pubs.acs.org/doi/pdf/10.1021/j100589a006
    # Car-Car can diverge from 1.34 to 1.443, Car-Nar from 1.325 to 1.39
    # C-Car : 1.51 to 1.53 very similar to C-C so is not counted
    # C-C : 1.53
    # C-N of lysine : 1.48
    # C-O for tyr is 1.36
    # C=O ar of ASP is from 1.24-1.29 --> use higher sigma
    # C-N of Arg (all 3) 1.32 to 1.35 fall within Car-Nar due to conjugation --> will be seperated based on one carbon with 3N

    bonds_dict = {'Car-Car': 1.40,
                  'Car-Nar': 1.34,
                  'C=O': 1.21,
                  'C-C': 1.54,
                  'C-N': 1.47,
                  'C-H': 1.09,
                  'C-O': 1.43,
                  'N-H': 0.99,
                  'C-S': 1.82,
                  'O-H': 0.98,
                  'S-H': 1.33,
                  'Car-O': 1.36
                  }

    counter_dict = {'Car-Car': 0,
                    'Car-Nar': 0,
                    'C=O': 0,
                    'C-C': 0,
                    'C-N': 0,
                    'C-H': 0,
                    'C-O': 0,
                    'N-H': 0,
                    'C-S': 0,
                    'O-H': 0,
                    'S-H': 0,
                    'Car-O': 0
                    }

    for idx, element in enumerate(inner_dist_matrix[row_idx]):
        if element > 1.6 and atoms_dict[idx] != 'S' and atoms_dict[row_idx] != 'S':
            break
        elif element > 1.9:
            #print('longer bond')
            break
            # if only python had a built in 'switch' functions this would be much neater...:
        if atoms_dict[row_idx] == 'H':
            return counter_dict  # not intersted in H for chain udentification skip them

        elif atoms_dict[row_idx] == 'N':
            # now find the index of element in the matrix before sorting to find which atom is it
            sorted_idx = sorted_idx_mat[row_idx][idx]
            if atoms_dict[sorted_idx] == 'H' and element < (bonds_dict['N-H'] + sigma_NH) and element > (
                    bonds_dict['N-H'] - sigma_NH):
                counter_dict['N-H'] += 1
                continue
            if atoms_dict[sorted_idx] == 'C' and element < (bonds_dict['Car-Nar'] + sigma) and element > (
                    bonds_dict['Car-Nar'] - sigma):
                counter_dict['Car-Nar'] += 1
                continue
            if atoms_dict[sorted_idx] == 'C' and element < (bonds_dict['C-N'] + sigma) and element > (
                    bonds_dict['C-N'] - sigma):
                counter_dict['C-N'] += 1
                continue

        elif atoms_dict[row_idx] == 'S':
            # now find the index of element in the matrix before sorting to find which atom is it
            sorted_idx = sorted_idx_mat[row_idx][idx]
            if atoms_dict[sorted_idx] == 'C' and element < (bonds_dict['C-S'] + sigma) and element > (
                    bonds_dict['C-S'] - sigma):
                counter_dict['C-S'] += 1
                continue
            if atoms_dict[sorted_idx] == 'H' and element < (bonds_dict['S-H'] + sigma) and element > (
                    bonds_dict['S-H'] - sigma):
                counter_dict['S-H'] += 1
                continue

        elif atoms_dict[row_idx] == 'O':
            # now find the index of element in the matrix before sorting to find which atom is it
            sorted_idx = sorted_idx_mat[row_idx][idx]
            if atoms_dict[sorted_idx] == 'H' and element < (bonds_dict['O-H'] + sigma) and element > (
                    bonds_dict['O-H'] - sigma):
                counter_dict['O-H'] += 1
                continue
            if atoms_dict[sorted_idx] == 'C' and element < (bonds_dict['Car-O'] + sigma) and element > (
                    bonds_dict['Car-O'] - sigma):
                counter_dict['Car-O'] += 1
                continue
            if atoms_dict[sorted_idx] == 'C' and element < (bonds_dict['C=O'] + sigma_CdO) and element > (
                    bonds_dict['C=O'] - sigma_CdO):
                counter_dict['C=O'] += 1
                continue
            if atoms_dict[sorted_idx] == 'C' and element < (bonds_dict['C-O'] + sigma) and element > (
                    bonds_dict['C-O'] - sigma):
                counter_dict['C-O'] += 1
                continue



        elif atoms_dict[row_idx] == 'C':
            sorted_idx = sorted_idx_mat[row_idx][idx]
            if atoms_dict[sorted_idx] == 'H' and element < (bonds_dict['C-H'] + sigma) and element > (
                    bonds_dict['C-H'] - sigma):
                counter_dict['C-H'] += 1
                continue
            if atoms_dict[sorted_idx] == 'N' and element < (bonds_dict['C-N'] + sigma) and element > (
                    bonds_dict['C-N'] - sigma):
                counter_dict['C-N'] += 1
                continue
            if atoms_dict[sorted_idx] == 'N' and element < (bonds_dict['Car-Nar'] + sigma) and element > (
                    bonds_dict['Car-Nar'] - sigma):
                counter_dict['Car-Nar'] += 1
                continue
            if atoms_dict[sorted_idx] == 'C' and element < (bonds_dict['C-C'] + sigma) and element > (
                    bonds_dict['C-C'] - sigma):
                counter_dict['C-C'] += 1
                continue
            if atoms_dict[sorted_idx] == 'C' and element < (bonds_dict['Car-Car'] + sigma) and element > (
                    bonds_dict['Car-Car'] - sigma):
                counter_dict['Car-Car'] += 1
                continue
            if atoms_dict[sorted_idx] == 'C' and element < (bonds_dict['C=O'] + sigma_CdO) and element > (
                    bonds_dict['C=O'] - sigma_CdO):
                counter_dict['C=O'] += 1
                continue
            if atoms_dict[sorted_idx] == 'C' and element < (bonds_dict['C-O'] + sigma) and element > (
                    bonds_dict['C-O'] - sigma):
                counter_dict['C-O'] += 1
                continue
            if atoms_dict[sorted_idx] == 'O' and element < (bonds_dict['Car-O'] + sigma) and element > (
                    bonds_dict['Car-O'] - sigma):
                counter_dict['Car-O'] += 1
            #     continue
            # if atoms_dict[sorted_idx] == 'S' and element < (bonds_dict['C-S'] + sigma) and element > (
            #         bonds_dict['C-S'] - sigma):
            #     counter_dict['C-S'] += 1
                continue

        else:  # if not Carbon or N or H for now:
            sys.exit('Input atom is not yet defined in this function, please write rivkakal@gmail.com')

    return counter_dict


def identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) -> str:
    '''identify residues by counting bonds types, returning residue name
        bond-lengths are teken from: 'Energy parameters in polypeptides. VII. Geometric parameters, partial atomic
        charges, nonbonded interactions ....etc' J.Phys.Chem. 1975,79,22,2361-2381 '''

    # create a dictinoary where matrix's row index correspond to the atom type!
    atoms_dict = {}
    keys = range(inner_dist_matrix.shape[0])
    values = []
    for idx, atom in enumerate(atoms_list):
        values.append(atom.get_id()[0])
    for i in keys:
        atoms_dict[i] = values[i]

    overall_counter_dict = {'Car-Car': 0,
                            'Car-Nar': 0,
                            'C-Car': 0,
                            'C=O': 0,
                            'C-C': 0,
                            'C-N': 0,
                            'C-H': 0,
                            'C-O': 0,
                            'N-H': 0,
                            'C-S': 0,
                            'O-H': 0,
                            'S-H': 0,
                            'Car-O': 0
                            }

    # now for every row in inner matrix, that is atom, assign bond types:
    for row in range(inner_dist_matrix.shape[0]):
        temp_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)
        overall_counter_dict = dict(Counter(overall_counter_dict) + Counter(temp_bonds_count_dict))
        # Check histidine with CA also
    # {'C-C': 4, 'C-H': 7, 'Car-Car': 2, 'Car-Nar': 8, 'N-H': 1} should be expected as C-H,N-H bonds are counted once,
    # but C-C,Car-Car,Car-Nar are counted twice!! so its equivalent ({'C-C': 2, 'C-H': 7, 'Car-Car': 1, 'Car-Nar': 4, 'N-H': 1} )
    his_dict = {'C-C': 2, 'C-H': 5, 'Car-Car': 2, 'Car-Nar': 8, 'N-H': 1}
    # for hsp there are two N-H bonds...
    res_flag = compare_dicts(overall_counter_dict, his_dict)

    if res_flag:
        return 'HIS'
    # now check positive his: HSP:
    hsp_dict = {'C-C': 2, 'C-H': 5, 'Car-Car': 2, 'Car-Nar': 8, 'N-H': 2}
    res_flag = compare_dicts(overall_counter_dict, hsp_dict)

    if res_flag:
        return 'HSP'

    phe_dict = {'C-C': 2, 'C-H': 8, 'Car-Car': 12}

    phe_flag = compare_dicts(overall_counter_dict, phe_dict)
    if phe_flag:
        return 'PHE'

    arg_dict = {'C-H': 3, 'C-N': 2, 'Car-Nar': 6, 'N-H': 5}
    arg_flag = compare_dicts(overall_counter_dict, arg_dict)

    if arg_flag:
        return 'ARG'


    lys_dict = {'C-N': 2, 'N-H': 3, 'C-H': 3}
    lys_flag = compare_dicts(overall_counter_dict, lys_dict)

    if lys_flag:
        return 'LYS'

    trp_dict = {'C-C': 2, 'Car-Car': 16, 'Car-Nar': 4, 'C-H': 8}
    trp_flag = compare_dicts(overall_counter_dict, trp_dict)
    if trp_flag:
        return 'TRP'

    # TODO check aromaticity in protein
    asn_dict = {'C-C': 2, 'C=O': 2, 'Car-Nar': 2, 'C-H': 3, 'N-H': 2}
    asn_flag = compare_dicts(overall_counter_dict, asn_dict)
    if asn_flag:
        return 'ASN'

    pro_dict = {'C-C': 6, 'C-N': 4, 'C-H': 8, 'N-H': 1}
    pro_flag = compare_dicts(overall_counter_dict, pro_dict)
    if pro_flag:
        return 'PRO'

    gly_dict = {'C-H': 4}  # in this analysis will also serve for all aliphatic only: ile, leu, val
    gly_flag = compare_dicts(overall_counter_dict, gly_dict)
    if gly_flag:
        return 'GLY'
    # TODO if ALA is before than it is recognized as ALA and not MET?!?!?
    met_dict = {'C-H': 6, 'C-S': 2}
    met_flag = compare_dicts(overall_counter_dict, met_dict)
    if met_flag:
        return 'MET'

    ala_dict = {'C-C': 2, 'C-H': 6}
    ala_flag = compare_dicts(overall_counter_dict, ala_dict)
    if ala_flag:
        return 'ALA'

    asp_dict = {'C-C': 2, 'C-H': 3, 'C=O': 4}
    # when non charged there is no 'O-H' ...
    asp_flag = compare_dicts(overall_counter_dict, asp_dict)
    if asp_flag:
        return 'ASP'

    cys_dict = {'C-S': 1, 'S-H': 1, 'C-H': 3}
    cys_flag = compare_dicts(overall_counter_dict, cys_dict)
    if cys_flag:
        return 'CYS'

    gln_dict = {'C-C': 2, 'C=O': 2, 'Car-Nar': 2, 'C-H': 3, 'N-H': 2}
    gln_flag = compare_dicts(overall_counter_dict, gln_dict)
    if gln_flag:
        return 'GLN'

    # glu is not modeled as aliphatic carbons here are not of interest (entropy is not measured here)
    # glu_dict = {'C-C': 2, 'C=O': 4, 'C-H': 3}
    # #'O-H' only when  non charged
    # glu_flag = compare_dicts(overall_counter_dict, glu_dict)
    # if glu_flag:
    #     return 'GLU'


    ser_dict = {'C-O': 2, 'O-H': 1, 'C-H': 3}
    ser_flag = compare_dicts(overall_counter_dict, ser_dict)
    if ser_flag:
        return 'SER'

    thr_dict = {'C-C': 2, 'C-O': 2, 'O-H': 1, 'C-H': 5}
    thr_flag = compare_dicts(overall_counter_dict, thr_dict)
    if thr_flag:
        return 'THR'

    tyr_dict = {'C-C': 2, 'Car-Car': 12, 'Car-O': 2, 'O-H': 1, 'C-H': 7}
    tyr_flag = compare_dicts(overall_counter_dict, tyr_dict)
    if tyr_flag:
        return 'TYR'

    else:
        raise ValueError(
            'Residue is not defined. Check for non necessary aliphatic carbons in tour structure (like CA...)')


def compare_dicts(current_dict: dict, to_compare_dict: dict) -> bool:
    # to compare dict has fewer values...
    for key in to_compare_dict:

        if key in current_dict:
            if not current_dict[key] == to_compare_dict[key]:
                res_flag = False
                break
            else:
                res_flag = True

    return res_flag


def assign_inner_matrix_bonds(inner_dist_matrix, atoms_list, sorted_idx_mat):
    'find covalent bonds according to residue type'
    # first create a dictinoary where matrix's row index correspond to the atom type!
    atoms_dict = {}
    keys = range(inner_dist_matrix.shape[0])
    values = []
    for idx, atom in enumerate(atoms_list):
        values.append(atom.get_id()[0])
    for i in keys:
        atoms_dict[i] = values[i]

    if identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) == 'HIS':  # could be HSD or HSE

        his_atoms_dict = {'CA': -1,
                          'CB': -1,
                          'CG': -1,
                          'ND1': -1,
                          'CE1': -1,
                          'NE2': -1,
                          'CD2': -1
                          }

        # CA_dict  = {'C-C': 1, 'C-H': 3} #here counting only once
        # CB_dict  = {'C-C': 2, 'C-H': 2}
        CB_dict = {'C-C': 1, 'C-H': 3}
        CG_dict = {'C-C': 1, 'Car-Nar': 1, 'Car-Car': 1}
        CE1_dict = {'Car-Nar': 2, 'C-H': 1}
        CD2_dict = {'Car-Car': 1, 'Car-Nar': 1, 'C-H': 1}
        # #now for HSD,HSE undefined H
        # NE2_dict = {'Car-Nar': 2}
        # ND1_dict = {'Car-Nar': 2}

        # now we need the row (index) of all atoms satisfying:
        for row in range(inner_dist_matrix.shape[0]):
            current_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)
            # if compare_dicts(current_bonds_count_dict, CA_dict):
            #         his_atoms_dict['CA'] = row
            #         continue
            if compare_dicts(current_bonds_count_dict, CB_dict):
                his_atoms_dict['CB'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CG_dict):
                his_atoms_dict['CG'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CE1_dict):
                his_atoms_dict['CE1'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CD2_dict):
                his_atoms_dict['CD2'] = row
                continue

        for row in range(inner_dist_matrix.shape[0]):
            if atoms_dict[row] == 'N':
                for idx, element in enumerate(inner_dist_matrix[row]):
                    if element > 1.6 and atoms_dict[idx] != 'S':
                        break

                    sorted_idx = sorted_idx_mat[row][idx]
                    if atoms_dict[sorted_idx] == 'C':
                        if sorted_idx == his_atoms_dict['CG']:  # this N is bonded to CG
                            his_atoms_dict['ND1'] = row
                        else:  # current N is not ND1 then
                            his_atoms_dict['NE2'] = row

        # find the HSD/HSE by which N has N-H (first)
        ND1_row_index = his_atoms_dict['ND1']
        if inner_dist_matrix[ND1_row_index][1] < 1.06 and inner_dist_matrix[ND1_row_index][
            1] > 0.95:  # this N has N-H bond
            his_type = 'HSD'
        else:
            his_type = 'HSE'
        return his_atoms_dict, his_type

    elif identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) == 'HSP':
        his_atoms_dict = {'CA': -1,
                          'CB': -1,
                          'CG': -1,
                          'ND1': -1,
                          'CE1': -1,
                          'NE2': -1,
                          'CD2': -1
                          }

        CB_dict = {'C-C': 1, 'C-H': 3}
        CG_dict = {'C-C': 1, 'Car-Nar': 1, 'Car-Car': 1}
        CE1_dict = {'Car-Nar': 2, 'C-H': 1}
        CD2_dict = {'Car-Car': 1, 'Car-Nar': 1, 'C-H': 1}

        # now we need the row (index) of all atoms satisfying:
        for row in range(inner_dist_matrix.shape[0]):
            current_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)

            if compare_dicts(current_bonds_count_dict, CB_dict):
                his_atoms_dict['CB'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CG_dict):
                his_atoms_dict['CG'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CE1_dict):
                his_atoms_dict['CE1'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CD2_dict):
                his_atoms_dict['CD2'] = row
                continue

        for row in range(inner_dist_matrix.shape[0]):
            if atoms_dict[row] == 'N':
                for idx, element in enumerate(inner_dist_matrix[row]):
                    if element > 1.6 and atoms_dict[idx] != 'S':
                        break

                    sorted_idx = sorted_idx_mat[row][idx]
                    if atoms_dict[sorted_idx] == 'C':
                        if sorted_idx == his_atoms_dict['CG']:  # this N is bonded to CG
                            his_atoms_dict['ND1'] = row
                        else:  # current N is not ND1 then
                            his_atoms_dict['NE2'] = row
        his_type = 'HSP'

        return his_atoms_dict, his_type

    elif identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) == 'PHE':
        phe_atoms_dict = {'CA': -1,
                          'CB': -1,
                          'CG': -1,
                          'CD1': -1,
                          'CE1': -1,
                          'CE2': -1,
                          'CD2': -1
                          }

        # CA_dict  = {'C-C': 1, 'C-H': 3} #here counting only once
        # CB_dict  = {'C-C': 2, 'C-H': 2}
        CB_dict = {'C-C': 1, 'C-H': 3}
        CG_dict = {'C-C': 1, 'Car-Car': 2}
        # CD1_dict = {'Car-Car':2, 'C-H': 1} #same for CD2,CE1,CE2 ...

        # now we need the row (index) of all atoms satisfying:
        for row in range(inner_dist_matrix.shape[0]):
            current_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)

            # if compare_dicts(current_bonds_count_dict, CA_dict):
            #         his_atoms_dict['CA'] = row
            #         continue
            if compare_dicts(current_bonds_count_dict, CB_dict):
                phe_atoms_dict['CB'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CG_dict):
                phe_atoms_dict['CG'] = row
                continue

        # look which carbons are bound to which...first find CD1,CD2 (doesn't matter really from symmetery issues) -->assign CE2,CD2 and then CZ
        # first time assign CD1 and CD2 . As we donot know if the CE2 attached to CD2 was defined before in the matrix,
        # we'll have to run over the same loop again as not to miss the assignent ot CE2 (by bonded to the now defined CD2)
        is_CD1_defined = False
        run1_finished = False
        run2_finished = False
        for row in range(inner_dist_matrix.shape[0]):
            if atoms_dict[row] == 'C' and not run1_finished:
                for idx, element in enumerate(inner_dist_matrix[row]):
                    if element > 1.6:
                        break
                    sorted_idx = sorted_idx_mat[row][idx]
                    if atoms_dict[sorted_idx] == 'C':
                        if sorted_idx == phe_atoms_dict['CG'] and not row == phe_atoms_dict['CB'] and not row == \
                                                                                                          phe_atoms_dict[
                                                                                                              'CG'] and not is_CD1_defined:  # this C is bonded to CG but is not CB (either CD1,CD2)
                            phe_atoms_dict['CD1'] = row
                            is_CD1_defined = True
                            continue
                        if sorted_idx == phe_atoms_dict['CG'] and not row == phe_atoms_dict[
                            'CB'] and is_CD1_defined and not row == phe_atoms_dict[
                            'CD1']:  # this C is bonded to CG but is not CB or CD1
                            phe_atoms_dict['CD2'] = row
                            run1_finished = True

        for row in range(inner_dist_matrix.shape[0]):
            if atoms_dict[row] == 'C' and not run2_finished:
                for idx, element in enumerate(inner_dist_matrix[row]):
                    if element > 1.6:
                        break
                    sorted_idx = sorted_idx_mat[row][idx]
                    if atoms_dict[sorted_idx] == 'C':
                        if sorted_idx == phe_atoms_dict['CD1'] and not row == phe_atoms_dict['CG'] and not row == \
                                                                                                           phe_atoms_dict[
                                                                                                               'CD1']:  # this C is bonded to CD1 but is not CG
                            phe_atoms_dict['CE1'] = row
                            continue
                        if sorted_idx == phe_atoms_dict['CD2'] and not row == phe_atoms_dict['CG'] and not row == \
                                                                                                           phe_atoms_dict[
                                                                                                               'CD2']:  # this C is bonded to CD1 but is not CG
                            phe_atoms_dict['CE2'] = row
                            run2_finished = True

        for row in range(inner_dist_matrix.shape[0]):
            if atoms_dict[row] == 'C':
                for idx, element in enumerate(inner_dist_matrix[row]):
                    if element > 1.6:
                        break
                    sorted_idx = sorted_idx_mat[row][idx]
                    if atoms_dict[sorted_idx] == 'C':
                        if sorted_idx == phe_atoms_dict['CE1'] and not row == phe_atoms_dict[
                            'CD1']:  # this C is bonded to CE1 but is not CD1
                            phe_atoms_dict['CZ'] = row
        res_type = 'PHE'
        return phe_atoms_dict, res_type

    elif identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) == 'LYS':
        lys_atoms_dict = {'CE': -1,
                          'NZ': -1
                          }

        CE_dict = {'C-N': 1, 'C-H': 3}
        NZ_dict = {'C-N': 1, 'N-H': 3}

        # now we need the row (index) of all atoms satisfying:
        for row in range(inner_dist_matrix.shape[0]):
            current_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)
            if compare_dicts(current_bonds_count_dict, CE_dict):
                lys_atoms_dict['CE'] = row
                continue
            if compare_dicts(current_bonds_count_dict, NZ_dict):
                lys_atoms_dict['NZ'] = row
                continue

        res_type = 'LYS'
        return lys_atoms_dict, res_type

    elif identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) == 'ARG':
        arg_atoms_dict = {'CD': -1,
                          'NE': -1,
                          'CZ': -1,
                          'NH2': -1,
                          'NH1': -1,
                          }

        CD_dict = {'C-N': 1, 'C-H': 3}
        NE_dict = {'Car-Nar': 1, 'C-N': 1, 'N-H': 1}
        CZ_dict = {'Car-Nar': 3}

        # NH2, NH1 are symmetrical ...

        # now we need the row (index) of all atoms satisfying:
        for row in range(inner_dist_matrix.shape[0]):
            current_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)
            if compare_dicts(current_bonds_count_dict, CD_dict):
                arg_atoms_dict['CD'] = row
                continue
            if compare_dicts(current_bonds_count_dict, NE_dict):
                arg_atoms_dict['NE'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CZ_dict):
                arg_atoms_dict['CZ'] = row
                continue

        # similar to PHE choose random
        is_NH1_defined = False
        run1_finished = False
        for row in range(inner_dist_matrix.shape[0]):
            if atoms_dict[row] == 'N' and not run1_finished:
                for idx, element in enumerate(inner_dist_matrix[row]):
                    if element > 1.6:
                        break
                    sorted_idx = sorted_idx_mat[row][idx]
                    if atoms_dict[sorted_idx] == 'C':
                        if sorted_idx == arg_atoms_dict['CZ'] and not row == arg_atoms_dict[
                            'NE'] and not is_NH1_defined:  # this N is bonded to CZ
                            arg_atoms_dict['NH1'] = row
                            is_NH1_defined = True
                            continue
                        if sorted_idx == arg_atoms_dict['CZ'] and not row == arg_atoms_dict[
                            'NE'] and is_NH1_defined and not row == arg_atoms_dict['NH1']:
                            arg_atoms_dict['NH2'] = row
                            run1_finished = True

        return arg_atoms_dict, 'ARG'

    elif identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) == 'TRP':
        trp_atoms_dict = {'CB': -1,
                          'CG': -1,
                          'CD1': -1,
                          'NE1': -1,
                          'CE2': -1,
                          'CD2': -1,
                          'CE3': -1,
                          'CZ3': -1,
                          'CH2': -1,
                          'CZ2': -1
                          }

        CB_dict = {'C-C': 1, 'C-H': 3}
        CG_dict = {'C-C': 1, 'Car-Car': 2}
        CD1_dict = {'Car-Car': 1, 'Car-Nar': 1, 'C-H': 1}
        NE1_dict = {'Car-Nar': 2, 'N-H': 1}
        CE2_dict = {'Car-Nar': 1, 'Car-Car': 2}
        CD2_dict = {'Car-Car': 3}
        # for CE3,CZ3,CH2,CZ2 same dictionaries, so following PHE identification
        # now we need the row (index) of all atoms satisfying:
        for row in range(inner_dist_matrix.shape[0]):
            current_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)
            if compare_dicts(current_bonds_count_dict, CB_dict):
                trp_atoms_dict['CB'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CG_dict):
                trp_atoms_dict['CG'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CD1_dict):
                trp_atoms_dict['CD1'] = row
                continue
            if compare_dicts(current_bonds_count_dict, NE1_dict):
                trp_atoms_dict['NE1'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CE2_dict):
                trp_atoms_dict['CE2'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CD2_dict):
                trp_atoms_dict['CD2'] = row
                continue

        for row in range(inner_dist_matrix.shape[0]):
            if atoms_dict[row] == 'C':
                for idx, element in enumerate(inner_dist_matrix[row]):
                    if element > 1.6:
                        break
                    sorted_idx = sorted_idx_mat[row][idx]
                    if atoms_dict[sorted_idx] == 'C':
                        if sorted_idx == trp_atoms_dict['CD2'] and not row == trp_atoms_dict['CD2'] and not row == \
                                                                                                            trp_atoms_dict[
                                                                                                                'CE2'] and not row == \
                                                                                                                               trp_atoms_dict[
                                                                                                                                   'CG']:  # this C is bonded to CD2 but is not CD2 and not CE2 and not CG
                            trp_atoms_dict['CE3'] = row
                            continue
                        if sorted_idx == trp_atoms_dict['CE2'] and not row == trp_atoms_dict['CE2'] and not row == \
                                                                                                            trp_atoms_dict[
                                                                                                                'CD2']:  # this C is bonded to CE2 but is not CD2 or CE2 and is not NE1 as it is carbon (not N)
                            trp_atoms_dict['CZ2'] = row
                            continue
        # second run after defining CZ1,CE3
        for row in range(inner_dist_matrix.shape[0]):
            if atoms_dict[row] == 'C':
                for idx, element in enumerate(inner_dist_matrix[row]):
                    if element > 1.6:
                        break
                    sorted_idx = sorted_idx_mat[row][idx]
                    if atoms_dict[sorted_idx] == 'C':
                        if sorted_idx == trp_atoms_dict['CE3'] and not row == trp_atoms_dict['CE3'] and not row == \
                                                                                                            trp_atoms_dict[
                                                                                                                'CD2']:
                            trp_atoms_dict['CZ3'] = row
                            continue
                        if sorted_idx == trp_atoms_dict['CZ2'] and not row == trp_atoms_dict['CZ2'] and not row == \
                                                                                                            trp_atoms_dict[
                                                                                                                'CE2']:
                            trp_atoms_dict['CH2'] = row
                            continue

        return trp_atoms_dict, 'TRP'

    elif identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) == 'ASN':
        asn_atoms_dict = {'CB': -1,
                          'CG': -1,
                          'OD1': -1,
                          'ND2': -1
                          }

        CB_dict = {'C-C': 1, 'C-H': 3}
        CG_dict = {'C=O': 1, 'Car-Nar': 1, 'C-C': 1}
        OD1_dict = {'C=O': 1}
        ND2_dict = {'Car-Nar': 2, 'N-H': 2}

        for row in range(inner_dist_matrix.shape[0]):
            current_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)
            if compare_dicts(current_bonds_count_dict, CB_dict):
                asn_atoms_dict['CB'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CG_dict):
                asn_atoms_dict['CG'] = row
                continue
            if compare_dicts(current_bonds_count_dict, OD1_dict):
                asn_atoms_dict['OD1'] = row
                continue
            if compare_dicts(current_bonds_count_dict, ND2_dict):
                asn_atoms_dict['ND2'] = row
                continue

        return asn_atoms_dict, 'ASN'

    # TODO issue with proline and the carbonyl missing imposing pseudo symmetry
    elif identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) == 'PRO':
        pro_atoms_dict = {'CA': -1,
                          'CB': -1,
                          'CG': -1,
                          'CD': -1,
                          'N': -1,
                          }

        N_dict = {'C-N': 2, 'N-H': 1}

        is_CA_defined = False
        first_run_finished = False
        # due to symmetry, analysis is similar to PHE
        for row in range(inner_dist_matrix.shape[0]):
            current_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)
            if compare_dicts(current_bonds_count_dict, N_dict):
                pro_atoms_dict['N'] = row
                break

        for row in range(inner_dist_matrix.shape[0]):
            if atoms_dict[row] == 'C' and not first_run_finished:
                for idx, element in enumerate(inner_dist_matrix[row]):
                    if element > 1.6:
                        break
                    sorted_idx = sorted_idx_mat[row][idx]
                    if atoms_dict[sorted_idx] == 'N':
                        if sorted_idx == pro_atoms_dict[
                            'N'] and not is_CA_defined:  # N bound carbon: assign first to be CA
                            pro_atoms_dict['CA'] = row
                            is_CA_defined = True
                            continue
                        if sorted_idx == pro_atoms_dict['N'] and is_CA_defined:  # N bound carbon: assign first to be CA
                            pro_atoms_dict['CD'] = row
                            first_run_finished = True
                            continue
        # second run after defining CA,CD
        for row in range(inner_dist_matrix.shape[0]):
            if atoms_dict[row] == 'C':
                for idx, element in enumerate(inner_dist_matrix[row]):
                    if element > 1.6:
                        break
                    sorted_idx = sorted_idx_mat[row][idx]
                    if atoms_dict[sorted_idx] == 'C':
                        if sorted_idx == pro_atoms_dict['CA'] and not row == pro_atoms_dict['CA']:
                            pro_atoms_dict['CB'] = row
                            continue
                        if sorted_idx == pro_atoms_dict['CD'] and not row == pro_atoms_dict['CD']:
                            pro_atoms_dict['CG'] = row
                            continue

        return pro_atoms_dict, 'PRO'

    elif identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) == 'GLY':
        gly_atoms_dict = {'CA': -1
                          }

        CA_dict = {'C-H': 4}

        for row in range(inner_dist_matrix.shape[0]):
            current_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)
            if compare_dicts(current_bonds_count_dict, CA_dict):
                gly_atoms_dict['CA'] = row
                break

        return gly_atoms_dict, 'GLY'

    elif identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) == 'ASP':
        asp_atoms_dict = {'CB': -1,
                          'CG': -1,
                          'OD1': -1,
                          'OD2': -1
                          }

        CB_dict = {'C-H': 3, 'C-C': 1}
        CG_dict = {'C-C': 1, 'C=O': 2}

        for row in range(inner_dist_matrix.shape[0]):
            current_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)
            if compare_dicts(current_bonds_count_dict, CB_dict):
                asp_atoms_dict['CB'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CG_dict):
                asp_atoms_dict['CG'] = row
                continue
        # symmetry in oxygens
        is_OD1_defined = False
        first_run_finished = False

        for row in range(inner_dist_matrix.shape[0]):
            if atoms_dict[row] == 'O' and not first_run_finished:
                for idx, element in enumerate(inner_dist_matrix[row]):
                    if element > 1.6:
                        break
                    sorted_idx = sorted_idx_mat[row][idx]
                    if atoms_dict[sorted_idx] == 'C':
                        if sorted_idx == asp_atoms_dict['CG'] and not is_OD1_defined:
                            asp_atoms_dict['OD1'] = row
                            is_OD1_defined = True
                            continue
                        if sorted_idx == asp_atoms_dict['CG'] and is_OD1_defined:
                            asp_atoms_dict['OD2'] = row
                            first_run_finished = True
                            continue

        return asp_atoms_dict, 'ASP'

    elif identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) == 'CYS':
        cys_atoms_dict = {'CB': -1,
                          'S': -1
                          }

        CB_dict = {'C-S': 1, 'C-H': 3}
        S_dict = {'C-S': 1, 'S-H': 1}

        for row in range(inner_dist_matrix.shape[0]):
            current_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)
            if compare_dicts(current_bonds_count_dict, CB_dict):
                cys_atoms_dict['CB'] = row
                continue
            if compare_dicts(current_bonds_count_dict, S_dict):
                cys_atoms_dict['S'] = row
                continue

        return cys_atoms_dict, 'CYS'

    elif identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) == 'MET':
        met_atoms_dict = {'CG': -1,
                          'SD': -1,
                          'CE': -1,
                          }

        S_dict = {'C-S': 2}

        for row in range(inner_dist_matrix.shape[0]):
            current_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)
            if compare_dicts(current_bonds_count_dict, S_dict):
                met_atoms_dict['SD'] = row
                break

        is_CG_defined = False
        first_run_finished = False

        for row in range(inner_dist_matrix.shape[0]):
            if atoms_dict[row] == 'C' and not first_run_finished:
                for idx, element in enumerate(inner_dist_matrix[row]):
                    if element < 1.6:  # note this is MET, looking for C-S bond > 1.8 length
                        continue
                    sorted_idx = sorted_idx_mat[row][idx]
                    if atoms_dict[sorted_idx] == 'S':
                        if sorted_idx == met_atoms_dict['SD'] and not is_CG_defined:
                            met_atoms_dict['CG'] = row
                            is_CG_defined = True
                            continue
                        if sorted_idx == met_atoms_dict['SD'] and is_CG_defined:
                            met_atoms_dict['CE'] = row
                            first_run_finished = True
                            continue

        return met_atoms_dict, 'MET'

    elif identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) == 'SER':
        ser_atoms_dict = {'CB': -1,
                          'OG': -1
                          }

        CB_dict = {'C-O': 1, 'C-H': 3}
        OG_dict = {'C-O': 1, 'O-H': 1}

        for row in range(inner_dist_matrix.shape[0]):
            current_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)
            if compare_dicts(current_bonds_count_dict, CB_dict):
                ser_atoms_dict['CB'] = row
                continue
            if compare_dicts(current_bonds_count_dict, OG_dict):
                ser_atoms_dict['OG'] = row
                continue

        return ser_atoms_dict, 'SER'

    elif identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) == 'THR':
        thr_atoms_dict = {'CB': -1,
                          'CG2': -1,
                          'OG1': -1
                          }

        CB_dict = {'C-O': 1, 'C-H': 2, 'C-C': 1}
        CG2_dict = {'C-C': 1, 'C-H': 3}
        OG1_dict = {'C-O': 1, 'O-H': 1}

        for row in range(inner_dist_matrix.shape[0]):
            current_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)
            if compare_dicts(current_bonds_count_dict, CB_dict):
                thr_atoms_dict['CB'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CG2_dict):
                thr_atoms_dict['CG2'] = row
                continue
            if compare_dicts(current_bonds_count_dict, OG1_dict):
                thr_atoms_dict['OG1'] = row
                continue

        return thr_atoms_dict, 'THR'

    # TODO HERE
    elif identify_sidechain(inner_dist_matrix, atoms_list, sorted_idx_mat) == 'TYR':
        tyr_atoms_dict = {'CB': -1,
                          'CG': -1,
                          'CD1': -1,
                          'CE1': -1,
                          'CZ': -1,
                          'OH': -1,
                          'CE2': -1,
                          'CD2': -1
                          }

        CB_dict = {'C-H': 3, 'C-C': 1}
        CG_dict = {'C-C': 1, 'Car-Car': 2}
        CZ_dict = {'Car-Car': 2, 'Car-O': 1}
        OH_dict = {'Car-O': 1, 'O-H': 1}

        for row in range(inner_dist_matrix.shape[0]):
            current_bonds_count_dict = find_bond_type(inner_dist_matrix, sorted_idx_mat, row, atoms_dict)
            if compare_dicts(current_bonds_count_dict, CB_dict):
                tyr_atoms_dict['CB'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CG_dict):
                tyr_atoms_dict['CG'] = row
                continue
            if compare_dicts(current_bonds_count_dict, CZ_dict):
                tyr_atoms_dict['CZ'] = row
                continue
            if compare_dicts(current_bonds_count_dict, OH_dict):
                tyr_atoms_dict['OH'] = row
                continue

        # symmetry in 4 ring atoms (in case of connectivity: 2 Car-Car and 1 C-H)
        is_CD1_defined = False
        run1_finished = False
        run2_finished = False
        for row in range(inner_dist_matrix.shape[0]):
            if atoms_dict[row] == 'C' and not run1_finished:
                for idx, element in enumerate(inner_dist_matrix[row]):
                    if element > 1.6:
                        break
                    sorted_idx = sorted_idx_mat[row][idx]
                    if atoms_dict[sorted_idx] == 'C':
                        if sorted_idx == tyr_atoms_dict['CG'] and not row == tyr_atoms_dict['CG'] and not row == \
                                                                                                          tyr_atoms_dict[
                                                                                                              'CB'] and not is_CD1_defined:  # this C is bonded to CG but is not CB (either CD1,CD2)
                            tyr_atoms_dict['CD1'] = row
                            is_CD1_defined = True
                            continue
                        if sorted_idx == tyr_atoms_dict['CG'] and not row == tyr_atoms_dict[
                            'CB'] and is_CD1_defined and not row == tyr_atoms_dict[
                            'CD1']:  # this C is bonded to CG but is not CB or CD1
                            tyr_atoms_dict['CD2'] = row
                            run1_finished = True

        for row in range(inner_dist_matrix.shape[0]):
            if atoms_dict[row] == 'C' and not run2_finished:
                for idx, element in enumerate(inner_dist_matrix[row]):
                    if element > 1.6:
                        break
                    sorted_idx = sorted_idx_mat[row][idx]
                    if atoms_dict[sorted_idx] == 'C':
                        if sorted_idx == tyr_atoms_dict['CD1'] and not row == tyr_atoms_dict['CG'] and not row == \
                                                                                                           tyr_atoms_dict[
                                                                                                               'CD1']:  # this C is bonded to CD1 but is not CG
                            tyr_atoms_dict['CE1'] = row
                            continue
                        if sorted_idx == tyr_atoms_dict['CD2'] and not row == tyr_atoms_dict['CG'] and not row == \
                                                                                                           tyr_atoms_dict[
                                                                                                               'CD2']:  # this C is bonded to CD1 but is not CG
                            tyr_atoms_dict['CE2'] = row
                            run2_finished = True

        return tyr_atoms_dict, 'TYR'



    # returns the atoms array and residue type for ligand change!!
    else:
        raise ValueError('residue is not recognized (check input structure bond lengths)')


######################################################### RENAMING PROCESS ##############################################################

def rename_current_atom(line_idx, line, new_atom, res_type):
    new_line = line
    if len(new_atom) == 3:
        curr_atom = new_line[12:15]
        new_line = new_line.replace(curr_atom, new_atom, 1)
    elif len(new_atom) == 2:
        curr_atom = new_line[12:15]
        new_line = new_line.replace(curr_atom, new_atom + ' ', 1)

    new_line = new_line.replace('LIG', res_type, 1)

    return new_line


def replace_numerical_with_greek(full_path: str, current_ligand_idx: str, previous_ligand_size: int, inner_dist_matrix,
                                 atoms_list, sorted_idx_mat) -> list:
    # for a given ligand id read text input and replace with the residue dictionary defenition
    lines_arr = []
    residue_dict, res_type = assign_inner_matrix_bonds(inner_dist_matrix, atoms_list, sorted_idx_mat)
    # print(residue_dict)
    # print(res_type)
    given_type = res_type
    key_list = list(residue_dict.keys())
    vall_list = list(residue_dict.values())

    if res_type == 'HSD' or res_type == 'HSE' or res_type == 'HSP':
        res_type = 'HIS'  # here don't care about his type, just to match with PDB ligand naming and work with already written functions

    with open(full_path, 'r') as rfile:
        for line_idx, line in enumerate(rfile):
            if line.startswith('HETATM'):
                splited_line = line.split()
                if splited_line[4] == current_ligand_idx:
                    if line[12] != ('H'):
                        position = vall_list.index(line_idx - previous_ligand_size)
                        new_atom = key_list[position]
                        new_line = rename_current_atom(line_idx - previous_ligand_size, line, new_atom, res_type)
                    else:
                        new_line = line
                        new_line = new_line.replace('LIG', res_type)
                    lines_arr.append(new_line)
    rfile.close()

    return lines_arr, given_type


############################################## GRREK OUTPUT SORTING ###################################################################
def create_greek_pdb_output(ligands_num: int, input_filename: str, dir_path, out_path):
    full_path = dir_path + input_filename
    final_arr = []
    added_str = ''  # this string will be added to file name
    geomtery = get_pdb_structure(full_path)
    geomtery_chain = create_chain_dict(geomtery)
    first_res = True
    previous_ligand_size = 0
    for i in range(1, ligands_num + 1):
        current_ligand_idx = i
        for chain_id, chain_obj in geomtery_chain.items():
            for residue in chain_obj.get_list():
                if residue.id[1] == current_ligand_idx:
                    atoms_list = residue.get_list()
                    size = get_ligand_size(full_path, current_ligand_idx)
                    #print(size)
                    mat = create_inner_dist_matrix(residue, size)
                    # mat = np.array(mat)
                    sorted_idx_mat = np.argsort(mat, axis=1)
                    mat.sort(axis=1)
                    temp_arr, res_type = replace_numerical_with_greek(full_path, str(current_ligand_idx),
                                                                      previous_ligand_size, mat, atoms_list,
                                                                      sorted_idx_mat)
                    previous_ligand_size = previous_ligand_size + size
                    final_arr.append(temp_arr)

            # if res_type is not None:
            #     added_str = added_str + res_type + '_' + str(i) + '_'

    # output_path = dir_path + 'residues' + added_str + '_' + input_filename[7:]
    output_path = out_path + 'trj_' + input_filename[7:]
    with open(output_path, 'w') as ofile:
        for ligand in final_arr:
            for line in ligand:
                ofile.write(line)
    ofile.close()
    return 'trj_' + input_filename[7:]
    # return final_arr


'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END of SECTION B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

# #call the output creating functions:
# #SECTION A:
# ligands_num = get_ligands_num(input_filename, dir_path)
# create_output(ligands_num, input_filename, dir_path)
# #SECTION B:
# final_input_filename = 'sorted_' + input_filename
# create_greek_pdb_output(ligands_num, final_input_filename, dir_path)
