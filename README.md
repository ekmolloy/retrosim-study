This repository contains data sets and scripts used in the study, ``Theoretical and practical considerations when using retroelement insertions to estimate species trees in the anomaly zone,'' available on [bioRxiv](https://doi.org/10.1101/2020.09.29.319038). The Python script to run ASTRAL_BP is [here](tools/run_astral_bp.py).

scripts directory
------------------
Contains bash scripts for 
+ [a](scripts/a_run_sp_methods.sh) : running parsimony methods
+ [b](scripts/b_run_astral_bp_methods.sh) : running ASTRAL_BP
+ [c](scripts/c_run_other_bp_methods.sh) : running ASTRID_BP and MDC_BP
+ [d](scripts/d_compute_species_tree_error.sh) : computing species tree error (writes CSV file)
+ [e](scripts/e_compute_species_tree_error_collapsed.sh) : collapsing branches based on some threshold for branch support and then computing species tree error (writes CSV file)
+ [f](scripts/f_grab_branch_lengths.sh) : extracting branch lengths and related data from newick strings (writes CSV file)
+ [g](scripts/g_run_astral_bp_for_rooted_Palaeognathae_trees.sh) : running ASTRAL_BP and re-rooting at chicken outgroup for Palaeognathae analysis
+ [h](scripts/h_analyze_Palaeognathae_5000_trees.sh)-[i](scripts/i_analyze_Palaeognathae_1000_trees.sh) : analyzing ASTRAL_BP trees from step (g) with 1000-5000 simulated retroelement insertions

Note that the scripts for the *ms* simulation and the species trees estimated using SDPquartets and Dollo parsimony are available on Dryad.

tools directory
----------------
Contains Python scripts called in bash scripts above

data directory
----------------
+ [astral-bp-trees](data/astral-bp-trees.tar.gz) : species trees estimated by running ASTRAL_BP
+ [astrid-bp-trees](data/astrid-bp-trees.tar.gz) : species trees estimated by running ASTRID_BP
+ [input-data-sets](data/input-data-sets.tar.gz) : retroelement insertion data sets simulated with ms
+ [mdc-bp-trees](data/mdc-bp-trees.tar.gz) : species trees estimated by running MDC_BP
+ [model-trees](data/model-trees.tar.gz) : true species trees used to simulate retroelement insertion data sets
+ [parsimony-trees](data/parsimony-trees.tar.gz) : species trees estimated by running variants of parsimony (unordered, Camin-Sokal, and Dollo)
+ [data-branch-info.csv](data/data-branch-info.csv) : branch lengths  estimated with ASTRAL_BP and related data (see [README](data/README_branch_info.md ))
+ [data-species-tree-error.csv](data/data-species-tree-error.csv) : species tree error for all methods and all data sets (see [README](data/README_species_tree_error.md))
+ [data-species-tree-error-collapsed.csv](data/data-species-tree-error-collapsed.csv) : species tree error for ASTRAL_BP collapsing branches based on their support values (see [README](data/README_species_tree_error.md))

summary directory
---------------------
Contains Python scripts for organizing data (see [csvs directory](summary/csvs) ), writing tables in Latex (see [tables directory](summary/tables)), and plotting data (see [plots directory](summary/plots))


external software
--------------------
Our experimental study uses several external packages:
+ ASTRAL v5.7.5 [[download]](https://github.com/smirarab/astral) [[paper]](https://doi.org/10.1186/s12859-018-2129-y)
+ ASTRID v2.2.1 [[download]](https://github.com/pranjalv123/ASTRID/releases/tag/2.2.1) [[paper]](https://doi.org/10.1186/1471-2164-16-S10-S3)
+ PAUP* v4a168 [[download]](http://paup.phylosolutions.com)
+ PhyloNet v3.8.2 [[download]](https://bioinfocs.rice.edu/phylonet/) [[paper]](https://doi.org/10.1093/sysbio/syy015)
+ PRANC [[download]](https://github.com/anastasiiakim/PRANC) [[paper]](https://doi.org/10.1093/molbev/msz305)
