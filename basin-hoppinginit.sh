#! /bin/bash

################################################################################################
#                                      Read the input.bh                                       #
################################################################################################

Simbolo_1=$(grep "cluster_ntyp" input.bh | cut -d "[" -f 2 | cut -d ":" -f 1)
Simbolo_2=$(grep "cluster_ntyp" input.bh | cut -d "[" -f 3 | cut -d ":" -f 1)
N_Simbolo_1=$(grep "cluster_ntyp" input.bh | cut -d "[" -f 2 | cut -d ":" -f 2 | cut -d "]" -f 1)
N_Simbolo_2=$(grep "cluster_ntyp" input.bh | cut -d "[" -f 3 | cut -d ":" -f 2 | cut -d "]" -f 1)
Nt1=$(grep "cluster_ntyp"  input.bh | cut -d " " -f 3 | cut -d "," -f 1 )
Nt2=$(grep "cluster_ntyp"  input.bh | cut -d " " -f 4)
XRange=$(grep "x_range"  input.bh | cut -d " " -f 3)
YRange=$(grep "y_range"  input.bh | cut -d " " -f 3)
ZRange=$(grep "z_range"  input.bh | cut -d " " -f 3)
ZVacuum=$(grep "z_vacuum"  input.bh | cut -d " " -f 3)
Pseudo_Dir=$( grep "pseudo_dir" input.bh | awk '{print $3}' )
pseudotype=$(grep "pseudo_type" input.bh | awk '{print $3}' ) 
step_width=$(grep "step_width" input.bh | awk '{print $3}')
Temperature=$(grep "temperature_K" input.bh | awk '{ print $3 }')
n=$(echo $Nt2  | wc -c  )
Nat=$(($N_Simbolo_1+$N_Simbolo_2)) #Numero de atomos del cluster
Ncore=$(grep "Ncore" input.bh | awk '{print $3}')
iteraciones=$(grep "iterations" input.bh | awk '{ print $3 }' )
path=$(grep "initialization_file" input.bh | awk '{ print $3 }' )
Npath=$(echo $path | wc -c )
m=1 #contador de coords rechazadas
Sel=$(grep "Selective"  input/POSCAR | wc -l )   #Determina si hay un selective dynamics

##############################################################################################

NPOSCAR=$(cat input/POSCAR | grep . | wc -l ) #Numero de lineas del poscar sin cluster
if [ $NPOSCAR -gt 5 ]
then
cd input
cat POSCAR | grep . >> aux
rm POSCAR
cat aux	>> POSCAR
rm aux
                                                       #-----------------------------------#
Species=$(head -6 POSCAR | tail -1)                    #                                   #
Numbers=$(head -7 POSCAR | tail -1)                    #                                   #
mv POSCAR aux                                          #            Estas lineas           #
head -5 aux >> POSCAR                                  #       concatenan los Símbolos     #
echo "$Species $Simbolo_1 $Simbolo_2" >> POSCAR        #            y números del          #
echo "$Numbers $N_Simbolo_1 $N_Simbolo_2 " >> POSCAR   #         clúster al POSCAR         #
tail -$(($NPOSCAR-7)) aux >> POSCAR                    #                                   #
rm aux                                                 #                                   #
                                                       #-----------------------------------#
else
cd input
cat POSCAR | grep . >> aux
rm POSCAR
cat aux >> POSCAR
rm aux

echo " $Simbolo_1 $Simbolo_2" >> POSCAR
echo "$N_Simbolo_1 $N_Simbolo_2 " >> POSCAR
echo "Cartesian">>POSCAR
fi

############################################################ PUEDE QUE ESTO YA NO SEA NECESARIO
#if [ $SelectiveDynamics = True ]
#then
#echo "Selective Dynamics" >> POSCAR
#echo "Cartesian">>POSCAR
#else
#echo "Cartesian">>POSCAR
#fi
#############################################################

################################################################################################
#                                   Generates POTCAR file                                      #
################################################################################################

head -6 POSCAR | tail -1 >> Species
tr -s '[:blank:]' '\n' < Species >> aux #Vuelve la fila una columna
grep . aux >> Symbols
rm aux
Ntyp=$(cat Symbols | wc -l )
rm Species
#Lee y guarda cada variable del simbolo
for((i=1;((i<$Ntyp+1));i++))
do
#S=$(eval "echo \$Simb_$i")
S=$(head -$i Symbols |tail -1 )
cat $Pseudo_Dir/$S$pseudotype/POTCAR>>POTCAR
done
#cd ..
rm Symbols


################################ GENERATES RUN.sh FILE #######################################

echo "mpirun -np $Ncore vasp_gam > output.out" >> run.sh
chmod +x run.sh

################################### Genera Matriz ############################################

head -5 POSCAR | tail -3 | awk '{print $1 " " $2 "  " $3}' >> ../programs/Matriz

##############################################################################################
cd ..

cp -r input $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2  ## OJO ACA: REVISAR COMO COPIA #Copia los archivos de input al directorio de trabajo
#python programs/RandomGenerator.py $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2/POSCAR   #Genera una estructura aleatoria y concatena a POSCAR

cd $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2
echo "  " >>aux
if [ $Npath -gt 8 ]
then   #Si existe path al archivo inicializador toma las coordenadas del mismo y las concatena al POSCAR
   tail -$Nat $path >>  POSCAR
else  #De otra forma  invoca al generador aleatorio
      if [ $n -gt 3 ] #Determina y corre de acuerdo con si es bimetálico  o monometálico #El 3 s arbitrario, solo determina si  [Au,3] es vacio o no
      then
         python ../programs/RandomGenerator.py aux $Nt1,$Nt2 $XRange $YRange $ZRange $ZVacuum
      else
         python ../programs/RandomGenerator.py aux $Nt1 $XRange $YRange $ZRange $ZVacuum
      fi
fi
for ((i=0;i<$Nat;i++))
do

echo " T  T  T">>din
done
tail -$Nat aux >> aux2
paste aux2 din >>  POSCAR

rm din aux aux2

#cd $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2  #De aqui en adelante trabaja en este directorio

mkdir rejected   #Crea el directorio de rechazados
./run.sh         #Crea el subproceso que corre VASP
contenido=$(grep "reached required accurcy" OUTCAR | wc -l )  # 0 si no converge, 1 si converge
#cd ..

#ACA INICIA LA  FASE 2 DE PRUEBA
while [[ $contenido != 1  ]];   #Mientras no converja la configuracion reinicia el calculo
do
echo " --> SCF failed. Starting again from randomly generated structure! "
      rm CHG CHGCAR DOSCAR EIGENVAL XDATCAR IBZKPT OSZICAR PCDAT REPORT WAVECAR *.xml CONTCAR POSCAR
      cp ../input/POSCAR POSCAR
echo "  " >>aux

     if [ $n -gt 3 ] #Determina y corre de acuerdo con si es bimetálico  o monometálico #El 3 s arbitrario, solo determina si  [Au,3] es vacio o no
      then
         python ../programs/RandomGenerator.py aux $Nt1,$Nt2 $XRange $YRange $ZRange $ZVacuum
      else
         python ../programs/RandomGenerator.py aux $Nt1 $XRange $YRange $ZRange $ZVacuum
      fi

for ((i=0;i<$Nat;i++))
do

echo " T  T  T">>din
done
tail -$Nat aux >> aux2
paste aux2 din >>  POSCAR
rm din aux aux2



./run.sh
contenido=$(grep "reached required accuracy" OUTCAR | wc -l )

#cd ..

done #Continua con el codigo si si convergio

#cd $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2
Energia=$(tail -1 OSZICAR | awk '{print $5 }')    #Extrae  la energia del OSZICA
echo "1 	$Energia " >> CONTCAR1
N=$(wc -l CONTCAR | awk '{ print $1 }' )   #Cuenta el numero de lineas que tiene CONTCAR
tail -$(($N-1)) CONTCAR >> CONTCAR1     #Mueve el CONTCAR a CONTCAR1
mv POSCAR POSCAR1    #Mueve POSCAR a POSCAR1
rm CHG CHGCAR DOSCAR EIGENVAL XDATCAR IBZKPT OSZICAR PCDAT REPORT WAVECAR *.xml CONTCAR

echo " "
echo " --> Relaxation of initial configuration: DONE! "
echo " "
echo "========================================================================================================="
echo "BH-DFT routine starts here! "
echo "Note: "
echo "For monometallic clusters: only random xyz moves will be applied "
echo "For bimetallic clusters  : 1 atomic swap will be performed after 10 moves "
echo "========================================================================================================="
echo " "

#COMIENZA FSE 3 E PRUEBA
i=2
while [ $i -lt $(($iteraciones+1)) ]  #Mientras i sea menor que el numero de iteraciones pedidas
do

let resto=$i%10   # Calcula el resto de i

   if [ $resto -eq 0 ]
   then #Aplica Swap

#*******************************************************************************************************************************************#

#      head -$(($NPOSCAR+3+$Nat)) CONTCAR$(($i-1)) >> aux2
 #     echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR
  #    head -$(($NPOSCAR+3)) aux2 | tail -$(($NPOSCAR+2)) >> POSCAR ##OJO ACA ESTO DARA PROLEMAS CON SELECTIVE DYNAMICS : PONER UN CONDICIONAL if [ $SelectiveDynamics = True ] ; then ; head -$NPOSCAR CONTCAR1 | tail -$(($NPOSCAR-1)) ; else head -8 CONTCAR1 | tail -7; fi
   #   tail -$Nat aux2 >> aux #Este contiene las coordenadas que leerá despues,
    #  rm aux2

#**************************************************************************************************#
 if [ $Sel -eq 1 ]
 then

    head -$(($NPOSCAR+$Nat)) CONTCAR$(($i+1)) >> aux2
    echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR
    head -$NPOSCAR CONTCAR1 | tail -$(($NPOSCAR-1)) >> POSCAR
    tail -$Nat aux2 >> aux #Este contiene las coordenadas que leerá despues,
    rm aux2

 else

    head -$(($NPOSCAR+3+$Nat)) CONTCAR$(($i-1)) >> aux2
    echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR
    head -8 CONTCAR1 | tail -7 >> POSCAR
    tail -$Nat aux2 >> aux #Este contiene las coordenadas que leerá despues,

 fi
#*******************************************************************************************************************************************#

tp=$(echo "$N_Simbolo_1
$N_Simbolo_2" | sort | head -1)


       for ((k=0;k<$tp;k++))
       do

      t=$(shuf -i 1-$N_Simbolo_1 -n 1)
      f=$(shuf -i 1-$N_Simbolo_2 -n 1) 
      m=$(($N_Simbolo_1+$f))
      lt=$(head -$t aux | tail -1 )
      lm=$(head -$m aux | tail -1 )

      for ((j=1;j<$(($Nat+1));j++))
      do
         if [ $j -eq $t ]
         then
            echo "$lm" >> POSCAR
         else
            if [ $j -eq $m ]
            then
               echo "$lt" >> POSCAR
            else
               head -$j aux | tail -1 >> POSCAR
            fi
          fi
       done
       rm aux

      done

#*******************************************************************************************************************************************#

   else #Aplica Move

#*******************************************************************************************************************************************#

#      head -$(($NPOSCAR+3+$Nat)) CONTCAR$(($i-1)) >> aux2
 #     echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR
  #    head -$(($NPOSCAR+3)) aux2 | tail -$(($NPOSCAR+2)) >> POSCAR  #OJO ACA: PROBLEMA CON SELECTIVE : IDEM
   #   tail -$Nat aux2 >> aux #Este contiene las coordenadas que leerá despues,
    #  rm aux2

#*******************************************************************************************************************************************#

 if [ $Sel -eq 1 ]
 then

    head -$(($NPOSCAR+$Nat)) CONTCAR$(($i+1)) >> aux2
    echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR
    head -$NPOSCAR CONTCAR1 | tail -$(($NPOSCAR-1)) >> POSCAR
    tail -$Nat aux2 >> aux #Este contiene las coordenadas que leerá despues,
    rm aux2

 else

    head -$(($NPOSCAR+3+$Nat)) CONTCAR$(($i-1)) >> aux2
    echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR
    head -8 CONTCAR1 | tail -7 >> POSCAR
    tail -$Nat aux2 >> aux #Este contiene las coordenadas que leerá despues,

 fi



#**********************************************************************************************************************#

      for ((j=1; j<$(($Nat+1)); j++))
      do

         dx=$(python ../programs/move.py $step_width)
         dy=$(python ../programs/move.py $step_width)
         dz=$(python ../programs/move.py $step_width)

         cd ../programs ; echo "$dx $dy $dz" | ./inverse > ../$Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2/kick 
         cd ../$Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2


         dxc=$(awk '{print $1}' kick)     #---------------------------------------------#
         dyc=$(awk '{print $2}' kick)     #    New (converted) coordinates of kick      #
         dzc=$(awk '{print $3}' kick)     #---------------------------------------------#

         x=$(head -$j aux | tail -1 | awk '{print $1}')
         y=$(head -$j aux | tail -1 | awk '{print $2}')
         z=$(head -$j aux | tail -1 | awk '{print $3}')

         echo "$(echo "$x+$xc" | bc ) $(echo "$y+$yc" | bc ) $(echo "$z+$dzc" | bc )" >> POSCAR

         # head -$j aux | tail -1 | awk '{print $1+$xc "  " $2+$yc "  " $3+$zc}' >> POSCAR

         rm kick

      done
      rm aux

   fi

   echo "Iteration $i of structure $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2"
   ./run.sh
   contenido=$(grep "reached required accuracy" OUTCAR | wc -l )

   while [[ $contenido != 1  ]];
   do
      echo " --> SCF failed. Starting again from randomly generated structure! "
      head -$(($NPOSCAR+3)) CONTCAR$(($i-1)) | tail -$(($NPOSCAR+2)) >> aux   #OJO: ESTODAR PROBLEMS CON Selective, agregar elmismo condicional

#*******************************************************************************************************************

 if [ $Sel -eq 1 ]
 then

    head -$(($NPOSCAR+$Nat)) CONTCAR$(($i+1)) >> aux2
    echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR
    head -$NPOSCAR CONTCAR1 | tail -$(($NPOSCAR-1)) >> POSCAR
    tail -$Nat aux2 >> aux #Este contiene las coordenadas que leerá despues,
    rm aux2

 else

    head -$(($NPOSCAR+3+$Nat)) CONTCAR$(($i-1)) >> aux2
    echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR
    head -8 CONTCAR1 | tail -7 >> POSCAR
    tail -$Nat aux2 >> aux #Este contiene las coordenadas que leerá despues,

 fi

#*******************************************************************************************************************
      # if [ $SelectiveDynamics = True ]
      # then
      # head -$NPOSCAR CONTCAR$(($i-1)) >> aux
      # else
      # head -$(($NPOSCAR+3)) CONTCAR$(($i-1)) >> aux
      # fi


      rm CHG CHGCAR DOSCAR EIGENVAL XDATCAR IBZKPT OSZICAR PCDAT REPORT WAVECAR *.xml CONTCAR POSCAR

      echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR

      cat aux >> POSCAR
#      python ../programs/replace_line.py POSCAR Direct Cartesian   ##OJO ACA: No reemplazar cart pordirect, usar inverse
      rm aux

      if [ $n -gt 3 ] #Determina y corre de acuerdo con si es bimetálico  o monometálico #El 3 s arbitrario, solo determina si  [Au,3] es vacio o no
      then
         python ../programs/RandomGenerator.py auxtoinvert $Nt1,$Nt2 $XRange $YRange $ZRange $ZVacuum
      else
         python ../programs/RandomGenerator.py auxtoinvert $Nt1 $XRange $YRange $ZRange $ZVacuum
      fi
echo "    " >> aux3
      for ((m=0;m<$Nat;m++))
      do
         head -$m auxtoinvert | tail -1 | ./../programs/inverse >> aux3
      done
for ((ij=0;ij<$Nat;ij++))
do

echo " T  T  T">>din
done
tail -$Nat aux3 >> aux4
paste aux4 din >>  POSCAR
rm din aux3 aux4


      rm auxtoinvert
      ./run.sh
      contenido=$(grep "reached required accuracy" OUTCAR | wc -l )

   done #Continua con el codigo si si convergio

   Energia=$(tail -1 OSZICAR | awk '{print $5 }')
   echo "$i	$Energia " >> CONTCAR$i
   N=$(wc -l CONTCAR | awk '{ print $1 }' )
   tail -$(($N-1)) CONTCAR >> CONTCAR$i
   mv POSCAR POSCAR$i
   rm CHG CHGCAR DOSCAR EIGENVAL XDATCAR IBZKPT OSZICAR PCDAT REPORT WAVECAR *.xml CONTCAR


   #COMIENZA Algoritmo de metropolis

   E0=$(head -1 CONTCAR$i | awk '{ print $1 }' )
   En=$(head -1 CONTCAR$(($i-1)) | awk '{ print $1 }' )

   accepted=$(python ../programs/metropolis.py $E0 $En $Temperature)

      if [ $accepted = true ]
      then
       	#Ha sido aceptado
	echo "--> Basin Hopping MC criteria: Energy accepted! "
	echo "--> Finished iteration $i"
	i=$(($i+1)) #Convergencia lograda, aumenta en 1 el contador
      else
	echo "--> Basin Hopping MC criteria: Energy rejected!"
	mv CONTCAR$i rejected/CONTCARrejected$m
        mv POSCAR$i rejected/POSCARrejected$m
        m=$(($m+1))
	let res i%10
#		if res=0
#			echo "--> Swap failed to converge"
#			swap_fail_flag=true
#		fi
	fi

   done

cd programs
rm Matriz
cd ..
cd input
rm run.sh POTCAR
cd ..
