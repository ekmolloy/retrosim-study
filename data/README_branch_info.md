CSV called *data-branch-info.csv* has the following columns:
---------------------------------------------------------------------


**DATA** : Retroelement insertion data set i.e.
+ 5taxa (see Figure 1A)
+ 6taxa (see Figure 1B
+ Palaeognathae (see Figure 1C)
+ 26taxa (see Figure 1D)
+ biological (Palaeognathaes retroelement insertion data set; see Table 1)


**REPL** : Replicate number
+ NA for biological data set
+ NA for model species trees (i.e. when MTHD = TRUE)
+ 1-25 for simulated data sets


**NRET** : Number of retroelement insertions in simulated data sets 
+ NA for biological data set
+ NA for model species trees (i.e. when MTHD = TRUE)
+ 100,000 for 5taxa and 6taxa simulated data sets
+ 10, 50, 100, 500, ... 100,000 for Palaeognathae and 26taxa simulated data sets



**MTHD** : 
+ TRUE
+ astral-bp = running ASTRAL_BP


**BIPA** : taxa on one side of branch (bipartition) -- used to identify the branch (note that side i.e. A vs B is chosen alphabetically)


**BIPB** : taxa on the other side of branch (bipartition)


**BRLN** : Branch length in coalescent units that is found after the semicolon in the newick string


**BRLN_MAP_GT** : Maximum a posteriori estimate of branch lengths, assuming gene trees as input (only available when running ASTRAL_BP with -t 2 option)


**BRLN_MLE_GT** : Maximum likelihood estimate of branch lengths, assuming gene trees as input (only available when running ASTRAL_BP with -t 2 option)


**BRLN_MLE_RI** : Maximum likelihood estimate of branch lengths, assuming retroelement insertions as input (only available when running ASTRAL_BP with -t 2 option)


**Q1** : normalized quartet support for dominant quartet topology (only available when running ASTRAL_BP with -t 2 option)


**Q2** : normalized quartet support for alternative quartet topology (only available when running ASTRAL_BP with -t 2 option)


**Q3** : normalized quartet support for other alternative quartet topology (only available when running ASTRAL_BP with -t 2 option)


**EN** : effective number of gene trees with information for resolving the branch, in this case retroelement insertions (only available when running ASTRAL_BP with -t 2 option)


**PP** : local posterior probability
