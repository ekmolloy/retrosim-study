#!/bin/bash

DATADIR="../data"
COMPARE="../tools/compare_two_trees.py"
CSV="../data/data-species-tree-error.csv"

if [ -e $CSV ]; then
    echo "$CSV already exists!"
    exit 1
fi

MODLS=( "5taxa" "6taxa" "26taxa" "Palaeognathae" )
REPLS=( $(seq 1 25) )
MTHDS=( "parsimony-strict" \
        "camin-sokal-strict" \
        "dollo-strict" \
        "astrid-bp" \
        "astral-bp" \
        "mdc-bp" \
        "mdc-bp-ur" )

echo "MODL,REPL,NRET,MTHD,NLEA,NINT_TRUE,NINT_ESTI,FN,FP,FNR,FPR,RF" > $CSV

for MODL in ${MODLS[@]}; do
    TRUE="$DATADIR/model-trees/${MODL}_model.tre"
    if [ $MODL == "5taxa" ] || [ $MODL == "6taxa" ]; then
        NRETS=( "100000" )
    else
        NRETS=( "10" "50" \
                "100" "500" \
                "1000" "5000" \
                "10000" "50000" \
                "100000" )
    fi
    for REPL in ${REPLS[@]}; do
        for NRET in ${NRETS[@]}; do
            for MTHD in ${MTHDS[@]}; do
                # Set path names
                if [ $MTHD == "astrid-bp" ] || [ $MTHD == "astral-bp" ]; then
                    INDIR="$DATADIR/${MTHD}-trees/$MODL"
                elif [ $MTHD == "mdc-bp" ] || [ $MTHD == "mdc-bp-ur" ]; then
                    INDIR="$DATADIR/mdc-bp-trees/$MODL"
                else
                    INDIR="$DATADIR/parsimony-trees/$MODL"
                fi

                if [ $NRET == "100000" ]; then
                        ESTI="$INDIR/$MTHD-Rep${REPL}.tre"
                else
                        ESTI="$INDIR/$MTHD-Rep${REPL}-${NRET}.tre"
                fi

                if [ -e $ESTI ]; then
                    if [ -z $(grep ";" $ESTI | awk '{print $1}') ]; then
                        DATA="NA,NA,NA,NA,NA,NA"
                    else
                        DATA=$(python $COMPARE -t $TRUE -e $ESTI)
                    fi
                else
                    DATA="NA,NA,NA,NA,NA,NA"
                fi

                echo "$MODL,$REPL,$NRET,$MTHD,$DATA" >> $CSV
            done
        done
    done
done
