#!/bin/bash
#BSUB -q q_residual
#BSUB -oo output
#BSUB -eo error
#BSUB -n 32
module load use.own
module load vasp/5.4.4
#BSUB -R "span[ptile=32]"
#BSUB -m "g3_a g3_b"
#./continue.sh 2> ERRORES
./pruebainit.sh 2> ERRORES > OUTPUT
#./Au1Ir1/aux.sh


