#!/bin/bash

INDIR="../data/input-data-sets"
OUTDIR="../data"

ASTRALSIESTA="../../software/ASTRAL-SIESTA/build/ASTRAL-SIESTA"
ASTRAL="../../software/ASTRAL-SIESTA/build/Astral/astral.5.7.5.jar"
ASTRID="../../software/ASTRID-osx"
PHYLONET="../../software/PhyloNet_3.8.2.jar"

CORRECT="python ../tools/correct_branch_lengths.py"
RUNMDCBP="python ../tools/run_mdc_bp.py"
RUNASTRALBP="python ../tools/run_astral_bp.py"


# Analyze simulated data sets with varying of numbers of retroelements
for MODL in "5taxa" "6taxa" "26taxa" "Palaeognathae"; do
    for REPL in `seq 1 25`; do
        # Analyze simulated data sets with 100,000 retroelements
        IN="$INDIR/$MODL/Rep${REPL}.bptrees"

        # Run ASTRAL-BP
        OUT="$OUTDIR/astral-bp-trees/$MODL/astral-bp-Rep${REPL}"
        if [ ! -e $OUT.log ]; then
            $RUNASTRALBP -ae $ASTRAL \
                         -i $INDIR/$MODL/Rep${REPL}.nex \
                         -gt -o $OUT &> $OUT.log
        fi

        # Run ASTRID
        OUT="$OUTDIR/astrid-bp-trees/$MODL/astrid-bp-Rep${REPL}"
        if [ ! -e $OUT.log ]; then
            $ASTRID -i $IN -o $OUT.tre &> $OUT.log
            rm $OUT.tre.1
        fi

        # Run MDC
        OUT="$OUTDIR/mdc-bp-trees/$MODL/mdc-bp-Rep${REPL}"
        if [ ! -e $OUT.tre ]; then
            echo "Estimating $OUT"
            $RUNMDCBP -e $PHYLONET -i $IN -o $OUT
        fi
        OUT="$OUTDIR/mdc-bp-trees/$MODL/mdc-bp-ur-Rep${REPL}"
        if [ ! -e $OUT.tre ]; then
            echo "Estimating $OUT"
                $RUNMDCBP -e $PHYLONET -u -i $IN -o $OUT
        fi

        if [ $MODL == "26taxa" ] || [ $MODL == "Palaeognathae" ]; then
            for NRET in "10" "50" "100" "500" "1000" "5000" "10000" "50000"; do
                IN="$INDIR/$MODL/Rep${REPL}-${NRET}.bptrees"
                head -n${NRET} $INDIR/$MODL/Rep${REPL}.bptrees > $IN

                # Run ASTRAL-BP
                OUT="$OUTDIR/astral-bp-trees/$MODL/astral-bp-Rep${REPL}-${NRET}"
                if [ ! -e $OUT.log ]; then
                    $RUNASTRALBP -ae $ASTRAL -n -i $IN -gt -o $OUT &> $OUT.log
                fi

                # Run ASTRID
                OUT="$OUTDIR/astrid-bp-trees/$MODL/astrid-bp-Rep${REPL}-${NRET}"
                if [ ! -e $OUT.log ]; then
                    $ASTRID -i $IN -o $OUT.tre &> $OUT.log
                    rm $OUT.tre.1    
                fi

                # Run MDC
                OUT="$OUTDIR/mdc-bp-trees/$MODL/mdc-bp-Rep${REPL}-${NRET}"
                if [ ! -e $OUT.tre ]; then
                    echo "Estimating $OUT"
                    $RUNMDCBP -e $PHYLONET -i $IN -o $OUT
                fi
                OUT="$OUTDIR/mdc-bp-trees/$MODL/mdc-bp-ur-Rep${REPL}-${NRET}"
                if [ ! -e $OUT.tre ]; then
                    echo "Estimating $OUT"
                    $RUNMDCBP -e $PHYLONET -u -i $IN -o $OUT
                fi

                rm $IN
            done
        fi
    done
done


# Analyze biological data sets
MODL="biological"
IN="$INDIR/$MODL/Palaeognathae.bptrees"
OUT="$OUTDIR/astral-bp-trees/$MODL/astral-bp"
if [ ! -e $OUT.log ]; then
    $RUNASTRALBP -ae $ASTRAL -n -i $IN -gt -o $OUT &> $OUT.log
fi
