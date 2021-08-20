#!/bin/bash

INDIR="../data/input-data-sets"
OUTDIR="../data"

PAUP="../../software/paup4a168_osx"
RUNPARSIMONY="python ../tools/run_parsimony.py"


# Analyze simulated data sets with varying of numbers of retroelements
for MODL in "5taxa" "6taxa" "26taxa" "Palaeognathae"; do
    OUTPATH="$OUTDIR/parsimony-trees/$MODL"
    if [ $MODL == "5taxa" ]; then
        OGS="Out"
        ARGS="-b"
    elif [ $MODL == "6taxa" ]; then
        OGS="Out"
        ARGS="-b"
    elif [ $MODL == "26taxa" ]; then
        OGS="TaxonU TaxonV TaxonW TaxonX TaxonY TaxonZ"
        ARGS=""
    else
        OGS="galGal"
        ARGS="-b"
    fi

    for REPL in `seq 1 25`; do
        IN="$INDIR/$MODL/Rep${REPL}.nex"

        # Run Unordered
        OUT="$OUTPATH/unordered-strict-Rep${REPL}.tre"
        if [ ! -e $OUT ]; then
            $RUNPARSIMONY -e $PAUP $ARGS -k 100000 -i $IN -g $OGS -o "$OUTPATH"
        fi

        # Run Dollo
        OUT="$OUTPATH/dollo-strict-Rep${REPL}.tre"
        if [ ! -e $OUT ]; then
            $RUNPARSIMONY -e $PAUP -p "dollo" $ARGS -k 100000 -i $IN -g $OGS -o "$OUTPATH"
        fi

        # Run Camin-Sokal
        OUT="$OUTPATH/camin-sokal-strict-Rep${REPL}.tre"
        if [ ! -e $OUT ]; then
            $RUNPARSIMONY -e $PAUP -p "camin-sokal" $ARGS -k 100000 -i $IN -g $OGS -o "$OUTPATH"
        fi

        # Varying number of RIs
        if [ $MODL == "26taxa" ] || [ $MODL == "Palaeognathae" ]; then
            for NRET in "10" "50" "100" "500" "1000" "5000" "10000" "50000"; do
                # Run Unordered
                OUT="$OUTPATH/unordered-strict-Rep${REPL}-$NRET.tre"
                if [ ! -e $OUT ]; then
                    $RUNPARSIMONY -e $PAUP $ARGS -k $NRET -i $IN -g $OGS -o "$OUTPATH"
                fi

                # Run Dollo
                OUT="$OUTPATH/dollo-strict-Rep${REPL}-$NRET.tre"
                if [ ! -e $OUT ]; then
                    $RUNPARSIMONY -e $PAUP -p "dollo" $ARGS -k $NRET -i $IN -g $OGS -o "$OUTPATH"
                fi

                # Run Camin-Sokal
                OUT="$OUTPATH/camin-sokal-strict-Rep${REPL}-$NRET.tre"
                if [ ! -e $OUT ]; then
                    $RUNPARSIMONY -e $PAUP -p "camin-sokal" $ARGS -k $NRET -i $IN -g $OGS -o "$OUTPATH"
                fi
            done
        fi
    done
done
