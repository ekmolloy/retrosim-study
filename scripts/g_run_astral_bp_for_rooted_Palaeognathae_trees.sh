#!/bin/bash

INDIR="../data/input-data-sets"
OUTDIR="../data"

ASTRALSIESTA="../../software/ASTRAL-SIESTA/build/ASTRAL-SIESTA"
ASTRAL="../../../software/ASTRAL-SIESTA/build/Astral/astral.5.7.5.jar"

RUNASTRALBP="python ../tools/run_astral_bp.py"


# Analyze simulated data sets with varying of numbers of retroelements
MODL="Palaeognathae"
for REPL in `seq 1 25`; do
    # Analyze simulated data sets with 100,000 retroelements
    IN="$INDIR/$MODL/Rep${REPL}.nex"

    # Run ASTRAL-BP
    #OUT="$OUTDIR/astral-bp-trees/$MODL/astral-bp-Rep${REPL}-rooted"
    #if [ ! -e $OUT.log ]; then
    #    $RUNASTRALBP -j $ASTRAL \
    #                 -g "galGal" \
    #                 -i $IN \
    #                 -o $OUT &> $OUT.log
    #fi

    #for NRET in "10" "50" "100" "500" "1000" "5000" "10000" "50000"; do
    for NRET in "1000" "5000"; do
        # Run ASTRAL-BP
        OUT="$OUTDIR/astral-bp-trees/$MODL/astral-bp-Rep${REPL}-${NRET}-rooted"
        if [ ! -e $OUT.log ]; then
            $RUNASTRALBP -j $ASTRAL \
                         -n $NRET \
                         -g "galGal" \
                         -i $IN \
                         -o $OUT \
                         &> $OUT.log
        fi
    done
done

# Analyze biological data sets
MODL="biological"
IN="$INDIR/$MODL/Palaeognathae.bptrees"
OUT="$OUTDIR/astral-bp-trees/$MODL/astral-bp-rooted"
if [ ! -e $OUT.log ]; then
    $RUNASTRALBP -newick -j $ASTRAL -i $IN -o $OUT -g "galGal" &> $OUT.log
fi
