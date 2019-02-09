Nat=$(($N_Simbolo_1+$N_Simbolo_2))
tail -$Nat inputcoordinatesfile >> aux
head -$(($(wc -l inputcoordinatesfile)-$Nat))>>POSCAR

#************************************************************************#
#             Ocupamos generar  numeros aleatorios con python            #
#************************************************************************#

for ((i=1; i<$(($Nat+1)); i++))
do

x=$(python move.py step_width)
y=$(python move.py step_width)
z=$(python move.py step_width)

head -$i aux | tail -1 | awk '{print $1+$x $2+$y $3+$z}' >> POSCAR

done
