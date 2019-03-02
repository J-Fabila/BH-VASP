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
NPOSCAR=$(wc -l input/POSCAR | awk '{ print $1}' ) #Numero de lineas del poscar sin cluster
Sel=$(grep "Selective"  input/POSCAR | wc -l )   #Determina si hay un selective dynamics


cd $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2  #Se mueve al directorio de trabajo
rm  CHG CHGCAR CONTCAR core* DOSCAR EIGENVAL IBZKPT OSZICAR OUTCAR PCDAT REPORT salida.out vasprun.xml WAVECAR XDATCAR 2> /dev/null  #Borra archivos residuales de la iteración  anterior
i=$(for i in $(ls CONTCAR* ); do head -1 $i | awk '{ print $1 }' ; done | sort -n  | tail -1 ) #Revisa en qué iteración se quedó
i=$(($i+1))                             #Actualiza el valor, inicia en la siguiente

cd rejected                #Ahora buscará en qué índice de rechazadas nos quedamos
m=$(ls CONTCAR* 2> /dev/null | wc -l  )
m=$(($m+1))

cd ..                                       #Regresamos al directorio de trabajo

while [ $i -lt $(($iteraciones+1)) ]           #Reanudamos el programa principal
do

let resto=$i%10

   if [ $Sel -eq 1 ]                 # Este if extrae las coordenadas de la iteracion anterior
   then                                          # dependiendo del caso: soportado o gas phase

      head -$(($NPOSCAR+$Nat)) CONTCAR$(($i-1)) >> aux2
      echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR
      head -$NPOSCAR CONTCAR1 | tail -$(($NPOSCAR-1)) >> POSCAR
      tail -$Nat aux2 >> aux                 #Este contiene las coordenadas que leerá despues
      rm aux2

   else

      head -$(($NPOSCAR+3+$Nat)) CONTCAR$(($i-1)) >> aux2
      echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR
      head -8 CONTCAR1 | tail -7 >> POSCAR
      tail -$Nat aux2 >> aux                 #Este contiene las coordenadas que leerá despues

   fi

   if [ $resto -eq 0 ]
   then #******************************** Aplica Swap ***************************************#


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
               echo "$lm" >> preposcar
            else
               if [ $j -eq $m ]
               then
                  echo "$lt" >> preposcar
               else
                  head -$j aux | tail -1 >> preposcar
               fi
            fi
         done
         rm aux

      done

   else #******************************** Aplica Move ***************************************#


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

         echo "$(echo "$x+$dxc" | bc ) $(echo "$y+$dyc" | bc ) $(echo "$z+$dzc" | bc )" >> preposcar

         rm kick

      done
      rm aux

   fi                                             # Cierra el move-swap. Continua el algoritmo

   if [ $Sel -eq 1 ] 
   then                                                   $

      for ((i=0;i<$Nat;i++))
      do

         echo " T  T  T">>din

      done


      tail -$Nat preposcar >> aux20
      paste aux20 din >>  POSCAR

      rm din aux20 preposcar
   else                      #De otra forma echa directamente las coordenadas de aux a POSCAR 

      tail -$Nat preposcar >> POSCAR
   fi


   echo "Iteration $i of structure $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2"

   ./run.sh
   contenido=$(grep "reached required accuracy" OUTCAR | wc -l )

   while [[ $contenido != 1  ]];
   do

      echo " --> SCF failed. Starting again from randomly generated structure! "

      if [ $Sel -eq 1 ]              # Este if extrae las coordenadas de la iteracion anterior
      then                                       # dependiendo del caso: soportado o gas phase

         head -$(($NPOSCAR+$Nat)) CONTCAR$(($i+1)) >> aux2
         echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR
         head -$NPOSCAR CONTCAR1 | tail -$(($NPOSCAR-1)) >> POSCAR
         tail -$Nat aux2 >> aux              #Este contiene las coordenadas que leerÃ¡ despues
         rm aux2

      else

         head -$(($NPOSCAR+3+$Nat)) CONTCAR$(($i-1)) >> aux2
         echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR
         head -8 CONTCAR1 | tail -7 >> POSCAR
         tail -$Nat aux2 >> aux              #Este contiene las coordenadas que leerÃ¡ despues

      fi

      rm CHG CHGCAR DOSCAR EIGENVAL XDATCAR IBZKPT OSZICAR PCDAT REPORT WAVECAR *.xml CONTCAR POSCAR

      echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR

      cat aux >> POSCAR
      rm aux

      echo "   " >> auxtoinvert            #Crea el archivo donde se guardarán las coordenadas

      if [ $n -gt 3 ] #Determina y corre de acuerdo con si es bimetÃ¡lico  o monometÃ¡lico
      then
         python ../programs/RandomGenerator.py auxtoinvert $Nt1,$Nt2 $XRange $YRange $ZRange $ZVacuum
      else
         python ../programs/RandomGenerator.py auxtoinvert $Nt1 $XRange $YRange $ZRange $ZVacuum
      fi


      for ((ma=0;ma<$Nat;ma++))
      do

         cd ../programs ;  head -$ma auxtoinvert | tail -1 | ./inverse >> ../$Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2/aux3 
         cd ../$Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2

      done

      if [ $Sel -eq 1 ]
      then                                                          #Si hay selective dynamics

         for ((i=0;i<$Nat;i++))
         do

            echo " T  T  T">>din

         done


      tail -$Nat aux3 >> aux4
      paste aux4 din >>  POSCAR

      rm din aux3 aux4

      else                    #De otra forma echa directamente las coordenadas de aux a POSCAR 

         tail -$Nat aux3 >> POSCAR
      fi

      rm auxtoinvert
      ./run.sh
      contenido=$(grep "reached required accuracy" OUTCAR | wc -l )

   done #Continua con el codigo si si convergio

   Energia=$(tail -1 OSZICAR | awk '{print $5 }')    #Extrae energia del OSZICAR
   echo "$i     $Energia " >> CONTCAR$i
   N=$(wc -l CONTCAR | awk '{ print $1 }' )
   tail -$(($N-1)) CONTCAR >> CONTCAR$i  #Renombra los archivos y les agrega la energia
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
#               if res=0
#                       echo "--> Swap failed to converge"
#                       swap_fail_flag=true
#               fi
        fi

   done

cd programs
rm Matriz
cd ..
cd input
rm run.sh POTCAR
cd ..
