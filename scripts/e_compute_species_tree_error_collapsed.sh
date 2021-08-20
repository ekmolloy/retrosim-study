#!/bin/bash

DATADIR="../data"
COMPARE="../tools/compare_two_trees.py"
CSV="../data/data-species-tree-error-collapsed.csv"

if [ -e $CSV ]; then
    echo "$CSV already exists!"
    exit 1
fi

MODLS=( "26taxa" "Palaeognathae" )
REPLS=( $(seq 1 25) )
NRETS=( "10" "50" "100" "500" "1000" \
        "5000" "10000" "50000" "100000" )
THRCS=( "0.4" "0.5" "0.6" \
        "0.7" "0.8" "0.9" )  # lowest support is 1/3 = 0.3333

echo "MODL,REPL,NRET,MTHD,THRC,NLEA,NINT_TRUE,NINT_ESTI,FN,FP,FNR,FPR,RF" > $CSV

for MODL in ${MODLS[@]}; do
    TRUE="$DATADIR/model-trees/${MODL}_model.tre"
    INDIR="$DATADIR/astral-bp-trees/$MODL"
    for REPL in ${REPLS[@]}; do
        for NRET in ${NRETS[@]}; do
            for THRC in ${THRCS[@]}; do
                if [ $NRET == "100000" ]; then
                        ESTI="$INDIR/astral-bp-Rep${REPL}.tre"
                else
                        ESTI="$INDIR/astral-bp-Rep${REPL}-${NRET}.tre"
                fi
                if [ -e $ESTI ]; then
                    if [ -z $(grep ";" $ESTI | awk '{print $1}') ]; then
                        DATA="NA,NA,NA,NA,NA,NA"
                    else
                        DATA=$(python $COMPARE -t $TRUE -e $ESTI -x $THRC)
                    fi
                else
                    DATA="NA,NA,NA,NA,NA,NA"
                fi
                echo "$MODL,$REPL,$NRET,astral-bp,$THRC,$DATA" >> $CSV
            done
        done
    done
done
