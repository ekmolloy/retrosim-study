CSV called *data-species-tree-error.csv* and *data-species-tree-error-collapsed.csv* has the following columns:
-------------------------------

**MODL** : Model species tree used to simulate retroelement insertion data set i.e.
+ 5taxa (see Figure 1A)
+ 6taxa (see Figure 1B
+ Palaeognathae (see Figure 1C)
+ 26taxa (see Figure 1D)


**REPL** : Replicate number
+ NA for biological data set
+ NA for model species trees (i.e. when MTHD = TRUE)
+ 1-25 for simulated data sets


**NRET** : Number of retroelement insertions in simulated data sets 
+ NA for biological data set
+ NA for model species trees (i.e. when MTHD = TRUE)
+ 100,000 for 5taxa and 6taxa simulated data sets
+ 10, 50, 100, 500, ... 100,000 for Palaeognathae and 26taxa simulated data sets


**MTHD** : Method used to estimate species tree from retroelement insertions
+ parsimony-strict   (data-species-tree-error.csv only)
+ camin-sokal-strict (data-species-tree-error.csv only)
+ dollo-strict       (data-species-tree-error.csv only)
+ astrid-bp          (data-species-tree-error.csv only)
+ astral-bp
+ mdc-bp             (data-species-tree-error.csv only)
+ mdc-bp-ur          (data-species-tree-error.csv only)


**THRC** (data-species-tree-error-collapsed.csv only) : Threshold for collapsing branches with low support
+ 0.4
+ 0.5
+ 0.6
+ 0.7
+ 0.8
+ 0.9


**NLEA** : Number of leaves in species tree


**NINT_TRUE** : Number of internal branches in unrooted true species tree (i.e. NLEA - 3)


**NINT_ESTI** : Number of internal branches in the unrooted estimated species tree


**FN** : False negatives i.e. internal branches in unrooted true species tree that are not in the unrooted estimated species tree


**FP** : False positives i.e. internal branches in the unrooted estimated species tree that are not in the unrooted true species tree


**FNR** : False negative rate i.e. FN / NINT_TRUE


**FPR** : False positive rate i.e. FP / NINT_TRUE


**RF** : (FN + FP) / (NINT_TRUE + NINT_TRUE)

