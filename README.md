# trypn
This repository contains python3 scripts for processing MS data from MASCOT-generated .csv files for determination of peptide, protein, and proteolysis parameters.

Each script has extended documentation for its use in its header.  These scripts perform the following general tasks:
    calc_cleavages.py
        determines N-terminal amino acid frequencies for cleavage events at both ends of assigned peptides
        
    comp_peptides_3.py
        determines Venn diagram parameters for overlapping peptides among three datasets
        
    comp_proteins_3.py
        determines Venn diagram parameters for overlapping proteins among three datasets
        
    compare_proteins_2.py
        determines Venn diagram parameters for overlapping proteins between two datasets
        
    compare_tryspin_trypn_peptides.py
        determines Venn diagram parameters for trypsin and Tryp-N digested peptides between two datasets
        
    compare_tryspin_trypn_peptides_no_mods.py
        determines Venn diagram parameters for trypsin and Tryp-N digested peptides between two datasets
        this script ignores modifications on peptides
