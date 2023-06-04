# pi_stacked_vs_cation_pi_QM_analysis
inclduing all files necessary to screen pi-pi and cation-pi pairs from a PDB dataset *.gz files. Also post-analysis of the calculated QM residue pairs parameters

# Pre analysis scripts
Searching all PDB high resolution datasets using for cation-pi pairs (HIS+,LYS,ARG) with (PHE,TYR,TRP,HIS0) using find_cationpi_pairs.py.
To cluster the pairs according to geometry use notebooks as **PHE_ARG_LYS_pre_clustering_analysis.ipynb** --> get a csv of all pairs to calculate at the QM level.
Make sure to include the following modules prior to the run **PlanesDistCalc.py, LoadStruct.py** (and modify the destination/input directories of dataset).
Similar search of pi-pi pairs according to their geometry is available using **specify_res_i_res_j_geometry_phi.py**, and for the clustering you may find **PHE_TYR_pre_analysis_clusteriny.ipynb**.

# QM calculations and Post-analysis notebooks
Representative structures are then cut into only relevant sidechain atoms (toluene for Phe, 4-methylimidazole for neutral HIS, etc' including only one additional carbon to the inteacting group)
Where HIS charge is involved, Hydrogens were added or removed using a PYMOL script, where coordinates are saved into xyz files.
QM calculations were running using ORCA 5.0.3 version's implemented DH: revDSD-PBE86-D4/QZ with CPCM(water) as solvent, up to 500 cycles of optimization.
obabel xyz to pdb was used initially to transform xyz to PDB format, yet residues were cut at sidechains atoms and thereford identified only as LIG instead of the correct residue name.
**LigandToResidue.py** is used to identify the correct LIG as residue. Is not fully implemented for other amino acids, though works for PHE,TYR,TRP,HIS(HSE,HSP),ARG and LYS.
Implementation of this module is available at **analyze_TYR_ARG.ipynb** or **analyze_PHE_TYR.py** along with other analysis including categorizing interactions.
