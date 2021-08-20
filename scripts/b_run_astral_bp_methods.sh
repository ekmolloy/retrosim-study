#!/bin/bash

INDIR="../data/input-data-sets"
OUTDIR="../data"

ASTRALSIESTA="../../software/ASTRAL-SIESTA/build/ASTRAL-SIESTA"
ASTRAL="../../software/ASTRAL-SIESTA/build/Astral/astral.5.7.5.jar"

RUNASTRALBP="python ../tools/run_astral_bp.py"


 Analyze simulated data sets with varying of numbers of retroelements
for MODL in "5taxa" "6taxa" "26taxa" "Palaeognathae"; do
    for REPL in `seq 1 25`; do
        # Analyze simulated data sets with 100,000 retroelements
        IN="$INDIR/$MODL/Rep${REPL}.nex"

        # Run ASTRAL-BP
        OUT="$OUTDIR/astral-bp-trees/$MODL/astral-bp-Rep${REPL}"
        if [ ! -e $OUT.log ]; then
            $RUNASTRALBP -j $ASTRAL \
                         -dosq \
                         -i $IN \
                         -o $OUT &> $OUT.log
        fi

        if [ $MODL == "26taxa" ] || [ $MODL == "Palaeognathae" ]; then
            for NRET in "10" "50" "100" "500" "1000" "5000" "10000" "50000"; do
                # Run ASTRAL-BP
                OUT="$OUTDIR/astral-bp-trees/$MODL/astral-bp-Rep${REPL}-${NRET}"
                if [ ! -e $OUT.log ]; then
                    $RUNASTRALBP -j $ASTRAL \
                                 -n $NRET \
                                 -dosq \
                                 -i $IN \
                                 -o $OUT \
                                 &> $OUT.log
                fi
            done
        fi
    done
done


# Analyze biological data sets
MODL="biological"
IN="$INDIR/$MODL/Palaeognathae.bptrees"
OUT="$OUTDIR/astral-bp-trees/$MODL/astral-bp"
if [ ! -e $OUT.log ]; then
    $RUNASTRALBP -newick -j $ASTRAL -i $IN -o $OUT &> $OUT.log
fi
