'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Details:
This Module contatins all functions related to creating a PDB structure object, checking chains and histidines a-priori!
Different functions are involved to handle non PDB information (casted quantum xyz->pdb characterized files)


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
from Bio.PDB import *

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PDB STRUCTURE LOADING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
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

def remove_irrelevant_chains(chain_dict: dict) -> dict:
    '''if any of the chains iterated is identical, remove the repeating units from calculations 
       to avoid biased data towards proteins with identical monomers (dimer, trimers...)
       Also checks wether the chain contains any HIS residue'''

    chains_to_remove = []
    for chain_id1, chain1 in chain_dict.items():
        hist_arr_1 = []
        for residue1 in chain1.get_list():
                    if residue1.resname == 'HIS':
                        hist_arr_1.append(residue1.get_id()[1])
        if len(hist_arr_1) < 1:
            chains_to_remove.append(chain_id1)
            continue
        for chain_id2, chain2 in chain_dict.items():
            if chain_id1 < chain_id2:
                hist_arr_2 = []            
                for residue2 in chain2.get_list():
                    if residue2.resname == 'HIS':
                        hist_arr_2.append(residue2.get_id()[1])
                if is_chain_identical(hist_arr_1, hist_arr_2):
                    chains_to_remove.append(chain_id2)


    for c2rm in chains_to_remove:
        print(c2rm)
        # TODO : not fully implemented yet
        chain_dict.pop(c2rm)
        
    return chain_dict


def check_histag(his_arr: list) -> list:
    """ check if there is a sequence of subsequent 6 histidines. If True, return their indices """
    if len(his_arr) > 6:
        for i, his in enumerate(his_arr):
            if i > len(his_arr) - 6:
                break
            histag_idxs = [his.get_id()]
            for j in range(1,7):
                if his.get_id()[1] + j == his_arr[j].get_id()[1]:
                    histag_idxs.append(his_arr[j].get_id()[1])
            if len(histag_idxs) > 5:
                return histag_idxs
    return []

def is_chain_identical(res_id_arr_a: list, res_id_arr_b: list) -> bool:
    """ checks if the given two residue occurances are the same size and their id's are identical """
    if len(res_id_arr_a) == len(res_id_arr_b):        
        if res_id_arr_a == res_id_arr_b:
            return True 
    return False


'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PDB STRUCTURE LOADING END OF SECTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

#############################################################################################################################################################
