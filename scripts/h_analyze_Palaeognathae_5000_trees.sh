#!/bin/bash

DATADIR="../data"
TRUE_TREE="$DATADIR/model-trees/Palaeognathae_model.tre"
INDIR="$DATADIR/astral-bp-trees/Palaeognathae"
NRET="5000"

COMPARE="../tools/compare_two_trees.py"
PREPARE2="../tools/prepare_to_plot_trees.py"
PREPARE4="../tools/prepare_for_az_test.py"
ISAZ="../tools/is_in_az.py"
FIND="../tools/find_incorrect_branch_with_max_support.py"

PHYLONET="../../../software/PhyloNet_3.8.2.jar"


#cat ../data/data-species-tree-error.csv | \
#    grep "Palaeognathae" | \
#    grep ",5000," | \
#    grep "astral-bp"

# Split into two groups based on correct and incorrect topologies
#Palaeognathae,3,5000,astral-bp,13,10,10,0,0,0.000000,0.000000,0.000000
#Palaeognathae,5,5000,astral-bp,13,10,10,0,0,0.000000,0.000000,0.000000
#Palaeognathae,6,5000,astral-bp,13,10,10,0,0,0.000000,0.000000,0.000000
#Palaeognathae,7,5000,astral-bp,13,10,10,0,0,0.000000,0.000000,0.000000
#Palaeognathae,9,5000,astral-bp,13,10,10,0,0,0.000000,0.000000,0.000000
#Palaeognathae,10,5000,astral-bp,13,10,10,0,0,0.000000,0.000000,0.000000
#Palaeognathae,11,5000,astral-bp,13,10,10,0,0,0.000000,0.000000,0.000000
#Palaeognathae,13,5000,astral-bp,13,10,10,0,0,0.000000,0.000000,0.000000
#Palaeognathae,14,5000,astral-bp,13,10,10,0,0,0.000000,0.000000,0.000000
#Palaeognathae,15,5000,astral-bp,13,10,10,0,0,0.000000,0.000000,0.000000
#Palaeognathae,16,5000,astral-bp,13,10,10,0,0,0.000000,0.000000,0.000000
#Palaeognathae,21,5000,astral-bp,13,10,10,0,0,0.000000,0.000000,0.000000
#Palaeognathae,23,5000,astral-bp,13,10,10,0,0,0.000000,0.000000,0.000000

#Palaeognathae,1,5000,astral-bp,13,10,10,1,1,0.100000,0.100000,0.100000
#Palaeognathae,2,5000,astral-bp,13,10,10,2,2,0.200000,0.200000,0.200000
#Palaeognathae,4,5000,astral-bp,13,10,10,1,1,0.100000,0.100000,0.100000
#Palaeognathae,8,5000,astral-bp,13,10,10,1,1,0.100000,0.100000,0.100000
#Palaeognathae,12,5000,astral-bp,13,10,10,1,1,0.100000,0.100000,0.100000
#Palaeognathae,17,5000,astral-bp,13,10,10,2,2,0.200000,0.200000,0.200000
#Palaeognathae,18,5000,astral-bp,13,10,10,1,1,0.100000,0.100000,0.100000
#Palaeognathae,19,5000,astral-bp,13,10,10,1,1,0.100000,0.100000,0.100000
#Palaeognathae,20,5000,astral-bp,13,10,10,1,1,0.100000,0.100000,0.100000
#Palaeognathae,22,5000,astral-bp,13,10,10,1,1,0.100000,0.100000,0.100000
#Palaeognathae,24,5000,astral-bp,13,10,10,1,1,0.100000,0.100000,0.100000
#Palaeognathae,25,5000,astral-bp,13,10,10,1,1,0.100000,0.100000,0.100000


ESTI_TREES=()
for REPL in `seq 1 25`; do
    ESTI_TREE="$INDIR/astral-bp-Rep${REPL}-${NRET}-rooted.tre"
    ESTI_TREES=( ${ESTI_TREES[@]} $ESTI_TREE )
done
python ../tools/compare_trees.py -t ${ESTI_TREES[@]}
#1 : 4 8 12 18 20 22 24
#2 : 
#3 : 5 6 7 9 10 11 13 14 15 16 21 23
#17 : 
#19 : 25


# NOW we prepare to plot some trees
mkdir -p $INDIR/trees-to-plot

# Case 1: correct topology
# + 13 replicates (3 5 6 7 9 10 11 13 14 15 16 21 23)
# Case 2: an incorrect topology
# + 8 replicates (1 4 8 12 18 20 22 24)
# Case 3: an incorrect topology
# + 2 replicates (19 and 25)
# Case 4: an incorrect topology
# + 1 replicate (2)
# Case 5: an incorrect topology
# + 1 replicate (17)

NAMEMAP="${INDIR}/name-map-to-plot.csv"

# Case 1 - 13 replicates (3 5 6 7 9 10 11 13 14 15 16 21 23)
INTREES=( "$INDIR/astral-bp-Rep3-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep5-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep6-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep7-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep9-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep10-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep11-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep13-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep14-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep15-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep16-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep21-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep23-${NRET}-rooted-t2.tre" )
OPREFIX="${INDIR}/trees-to-plot/astral-bp-Rep3-${NRET}-rooted-to-plot"
if [ ! -e $OPREFIX.tre ]; then
    python $PREPARE2 \
        -t ${INTREES[@]} \
        -m $NAMEMAP \
        -o $OPREFIX
fi

# Case 2 - 8 replicates (1 4 8 12 18 20 22 24)
INTREES=( "$INDIR/astral-bp-Rep1-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep4-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep8-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep12-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep18-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep20-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep22-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep24-${NRET}-rooted-t2.tre" )
OPREFIX="${INDIR}/trees-to-plot/astral-bp-Rep1-${NRET}-rooted-to-plot"
if [ ! -e $OPREFIX.tre ]; then
    python $PREPARE2 \
        -t ${INTREES[@]} \
        -m $NAMEMAP \
        -o $OPREFIX
fi

# Case 3
INTREES=( "$INDIR/astral-bp-Rep19-${NRET}-rooted-t2.tre" \
          "$INDIR/astral-bp-Rep25-${NRET}-rooted-t2.tre" )
OPREFIX="${INDIR}/trees-to-plot/astral-bp-Rep19-${NRET}-rooted-to-plot"
if [ ! -e $OPREFIX.tre ]; then
    python $PREPARE2 \
        -t ${INTREES[@]} \
        -m $NAMEMAP \
        -o $OPREFIX
fi

# Case 4
INTREES=( "$INDIR/astral-bp-Rep2-${NRET}-rooted-t2.tre"  )
OPREFIX="${INDIR}/trees-to-plot/astral-bp-Rep2-${NRET}-rooted-to-plot"
if [ ! -e $OPREFIX.tre ]; then
    python $PREPARE2 \
        -t ${INTREES[@]} \
        -m $NAMEMAP \
        -o $OPREFIX
fi

# Case 5
INTREES=( "$INDIR/astral-bp-Rep17-${NRET}-rooted-t2.tre"  )
OPREFIX="${INDIR}/trees-to-plot/astral-bp-Rep17-${NRET}-rooted-to-plot"
if [ ! -e $OPREFIX.tre ]; then
    python $PREPARE2 \
        -t ${INTREES[@]} \
        -m $NAMEMAP \
        -o $OPREFIX
fi


# NOW we check which trees are in the anomaly zone...
mkdir -p $INDIR/trees-for-az-test

NAMEMAP="${INDIR}/name-map-for-az-test.csv"

# The true tree is
# (((((A,B),C),D),ostrich),chicken)
#
# where
#
# A  = kiwi
#    = ((little_spotted_kiwi,great_spotted_kiwi),Okarito_brown_kiwi)
# B  = emu+sc
#    = (emu,southern_cassowary)
# AB = (A,B)
#    = (((little_spotted_kiwi,great_spotted_kiwi),Okarito_brown_kiwi),(emu,southern_cassowary))
# C  = rhea
#    = (lesser_rhea,greater_rhea)
# D  = tinamou
#    = ((thicket_tinamou,white_throated_tinamou),(Chilean_tinamou,elegant_crested_tinamou))

# Then looking at the 5 topologies recovered

# Case 1: 
# st = gt1 = (((A,B),C),D)
# gt = gt2

# Case 2:
# st = gt2 = ((A,B),(C,D))
# gt = probably not in AZ bc symmetric

# Case 3:
# st = gt3 = (((A,B),D),C)
# gt = gt2

# Case 4: 
# st = gt4 = (((A,D),C),B)
# gt = gt6 = ((A,D),(B,C))

# Case 5:
# st = gt5 = (((A,D),B),C)
# gt = gt6


# Prepare species tree (and first gene tree) for AZ test...
GTREES=()
for i in `seq 1 6`; do
  GTREES=( ${GTREES[@]} "${INDIR}/trees-for-az-test/gt${i}.tre" \ )
done

for REPL in `seq 1 25`; do
    STREE="$INDIR/astral-bp-Rep${REPL}-${NRET}-rooted.tre"
    OUTNEX="${INDIR}/trees-for-az-test/astral-bp-Rep${REPL}-${NRET}-rooted.nex"
    OUTLOG="${INDIR}/trees-for-az-test/astral-bp-Rep${REPL}-${NRET}-rooted.log"
    OUTTXT="${INDIR}/trees-for-az-test/astral-bp-Rep${REPL}-${NRET}-rooted.out"

    if [ ! -e $OUTNEX ]; then
        python $PREPARE4 \
            -s $STREE \
            -g ${GTREES[@]} \
            -m $NAMEMAP \
            -o $OUTNEX
    fi

    if [ ! -e $OUTLOG ]; then
        java -jar $PHYLONET $OUTNEX &> $OUTLOG
    fi

    if [ ! -e $OUTTXT ]; then
        cat $OUTLOG | python $ISAZ &> $OUTTXT
    fi
done


# Check for some other things...
python $FIND -t $TRUE_TREE -e ${ESTI_TREES[@]}
#Avg 0.557099
#Std 0.128334
#Min 0.370668
#Max 0.827469

