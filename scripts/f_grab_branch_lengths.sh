#!/bin/bash

DATADIR="../data"
GRAB="../tools/grab_branch_info.py"

CSV="$DATADIR/data-branch-info.csv"
if [ -e $CSV ]; then
	echo "$CSV already exists!"
	exit
fi
echo "DATA,REPL,NRET,MTHD,BIPA,BIPB,BRLN,BRLN_MAP_GT,BRLN_MLE_GT,BRLN_MLE_RI,Q1,Q2,Q3,EN,PP" > $CSV

python $GRAB \
    -i $DATADIR/astral-bp-trees/biological/astral-bp-t2.tre \
    -p "biological,NA,NA,astral-bp" >> $CSV

NRETS=( "10" "50" "100" "500" \
        "1000" "5000" "10000" "50000" )


for MODL in "5taxa" "6taxa" "26taxa" "Palaeognathae"; do
    python $GRAB \
        -i $DATADIR/model-trees/${MODL}_model.tre \
        -p "$MODL,NA,NA,true" >> $CSV

    for REPL in `seq 1 25`; do
        python $GRAB \
                -i $DATADIR/astral-bp-trees/$MODL/astral-bp-Rep${REPL}-t2.tre \
                -p "$MODL,$REPL,100000,astral-bp" >> $CSV

        if [ $MODL == "26taxa" ] || [ $MODL == "Palaeognathae" ]; then
            for NRET in ${NRETS[@]}; do
                python $GRAB \
                    -i $DATADIR/astral-bp-trees/$MODL/astral-bp-Rep${REPL}-${NRET}-t2.tre \
                    -p "$MODL,$REPL,$NRET,astral-bp" >> $CSV
            done
        fi
    done
done
