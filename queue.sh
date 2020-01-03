#!/bin/bash
#BSUB -q q_residual
#BSUB -oo output
#BSUB -eo error
#BSUB -n 128
module load use.own
module load vasp/5.4.4
#BSUB -R "span[ptile=32]"
#BSUB -m "g3_a g3_b"
#BSUB -J Ir5

#./continue.sh 2>> Ir5.err >> Ir5.out
./basin_hopping.sh 2>> Ir5.err >> Ir5.out


