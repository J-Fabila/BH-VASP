# BH-VASP

BH-VASP
│   input.bh
│   queue.sh
|   basin_hopping.sh
|   continue.sh
└───input
│        INCAR
│        KPOINTS
|        POSCAR
│        to_xyz
│        to_vasp
└───programs
         atomicpp.h
         inverse
         kickpp
         metropolis.py
         move.py
         RandomGenerator
         SRandomGenerator

file_name (output directory)
│   CONTCAR1, CONTCAR2, ... , CONTCARn
|   POSCAR1, POSCAR2, ..., POSCARn
|   energies.txt
|   sorted.txt
|   summary.txt
|   run.sh
|   INCAR, POTCAR, KPOINTS
|   CHG, CHGCAR, EIGENVAL, ... , IBZKPT
└───input (copy of the original)
└───rejected (by the MC criterion)
         CONTCARrejected12
         POSCARrejected12
         CONTCARrejected23
         POSCARrejected23
         
         
         
