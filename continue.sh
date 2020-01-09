Simbolo_1=$(grep "cluster_ntyp" input.bh | cut -d "[" -f 2 | cut -d ":" -f 1)
Simbolo_2=$(grep "cluster_ntyp" input.bh | cut -d "[" -f 3 | cut -d ":" -f 1) 2>/dev/null
N_Simbolo_1=$(grep "cluster_ntyp" input.bh | cut -d "[" -f 2 | cut -d ":" -f 2 | cut -d "]" -f 1)
N_Simbolo_2=$(grep "cluster_ntyp" input.bh | cut -d "[" -f 3 | cut -d ":" -f 2 | cut -d "]" -f 1) 2>/dev/null
Nt1=$(grep "cluster_ntyp"  input.bh | cut -d " " -f 3 | cut -d "," -f 1 )
Nt2=$(grep "cluster_ntyp"  input.bh | cut -d " " -f 4) 2>/dev/null
n=$(grep "cluster_ntyp"  input.bh | awk '{print $4}' | wc -c )  #Este es un criterio para determinar si es bimetálico o no

randomness=$(grep "randomness" input.bh | awk '{print $3}')
kick=$(grep "kick_type" input.bh | awk '{print $3}')
file_name=$(grep "file_name" input.bh | awk '{print $3}')

Pseudo_Dir=$( grep "pseudo_dir" input.bh | awk '{print $3}' )
pseudotype=$(grep "pseudo_type" input.bh | awk '{print $3}' )
step_width=$(grep "step_width" input.bh | awk '{print $3}')
Temperature=$(grep "temperature_K" input.bh | awk '{ print $3 }')
if [ $n -gt 3 ]
then
Nat=$(($N_Simbolo_1+$N_Simbolo_2)) #Numero de atomos del cluster
else
Nat=$(echo $N_Simbolo_1)
fi
#Ncore=$(grep "Ncore" input.bh | awk '{print $3}')
Ncore=$(grep " -n " queue.sh | awk '{print $3}' )
iteraciones=$(grep "iterations" input.bh | awk '{ print $3 }' )
path=$(grep "initialization_file" input.bh | awk '{ print $3 }' )
Npath=$(echo $path | wc -c )
m=1 #contador de coords rechazadas
Sel=$(grep "Selective"  input/POSCAR | wc -l )   #Determina si hay un selective dynamics
NPOSCAR=$(cat input/POSCAR | grep . | wc -l ) #Numero de lineas del poscar sin cluster
swap_step=10      # SWAP STEP
cd $file_name
rm POSCAR CHG CHGCAR CONTCAR core* DOSCAR EIGENVAL IBZKPT OSZICAR OUTCAR PCDAT REPORT salida.out vasprun.xml WAVECAR XDATCAR 2> /dev/null
i=$(for i in $(ls CONTCAR* ); do head -1 $i | awk '{ print $1 }' ; done | sort -n  | tail -1 )
i=$(($i+1))
cd rejected

m=$(ls CONTCAR* 2> /dev/null | wc -l  )
m=$(($m+1))
cd ../..

################################################################################

if [ $NPOSCAR -gt 8 ]                        #Determina si el sistema es gas phase o soportado
then                   #Si es soportado agrega los simbolos atomicos del cluster y sus numeros
   #Si es soportado además lee el cristal del POSCAR y establece los xrange,yrange y zrange.
   cd input
   cp ../programs/read_crystal .
   ./read_crystal > rangos
   rm read_crystal
   xmin=$(awk '{print $1}' rangos)
   xmax=$(awk '{print $2}' rangos)
   ymin=$(awk '{print $3}' rangos)
   ymax=$(awk '{print $4}' rangos)
   zmin=$(awk '{print $5}' rangos)
   zmax=$(awk '{print $6}' rangos)
   z_vac_min=$(awk '{print $7}' rangos)
   z_vac_max=$(awk '{print $8}' rangos)
   rm rangos
   cd ..
else          #Si es gas phase prepara el POSCAR, agrega los simbolos atomicos y el "Cartesian"
   cd programs
   ./read_matriz Matriz_$file_name > rangos
   xmin=$(awk '{print $1}' rangos)
   xmax=$(awk '{print $2}' rangos)
   ymin=$(awk '{print $3}' rangos)
   ymax=$(awk '{print $4}' rangos)
   zmin=$(awk '{print $5}' rangos)
   zmax=$(awk '{print $6}' rangos)
   z_vac_min=$(awk '{print $7}' rangos)
   z_vac_max=$(awk '{print $8}' rangos)
   rm rangos
   cd ..
fi

cd $file_name                                           #Nos mueve a ese directorio de trabajo#

################################################################################

while [ $(($i+$m)) -lt $(($iteraciones+1)) ]
do

   let resto=$i%$swap_step                  # Calcula el resto de i  para hacer los pasos swap

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
      rm aux2

   fi


   if [ $resto -eq 0 ]
   then #******************************** Aplica Swap ***************************************#

if [ $n -gt 3 ]
then
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
         mv preposcar aux
      done
      mv aux preposcar
else
   mv aux preposcar
fi


   else #******************************** Aplica Move ***************************************#
      if [ $kick -eq 0 ]
      then
         for ((j=1; j<$(($Nat+1)); j++))
         do

         dx=$(python ../programs/move.py $step_width)
         dy=$(python ../programs/move.py $step_width)
         dz=$(python ../programs/move.py $step_width)

         cd ../programs ; echo "$dx $dy $dz" | ./inverse Matriz_$file_name > ../$file_name/kick
echo "=====================kick========================="
cat ../$file_name/kick
echo "=================================================="
         cd ../$file_name

         dxc=$(awk '{print $1}' kick)     #---------------------------------------------#
         dyc=$(awk '{print $2}' kick)     #    New (converted) coordinates of kick      #
         dzc=$(awk '{print $3}' kick)     #---------------------------------------------#

         x=$(head -$j aux | tail -1 | awk '{print $1}')
         y=$(head -$j aux | tail -1 | awk '{print $2}')
         z=$(head -$j aux | tail -1 | awk '{print $3}')

         echo "$(echo "$x+$dxc" | bc ) $(echo "$y+$dyc" | bc ) $(echo "$z+$dzc" | bc )" >> preposcar

echo "============== OPERACIONES ================="
echo "$x+$dxc=$(echo "$x+$dxc" | bc )    $y+$dyc= $(echo "$y+$dyc" | bc )   $z+$dzc=$(echo "$z+$dzc" | bc )"
echo "============================================"
cat preposcar
         rm kick
echo "MOVE Performed"
         done
      else
         #Aplica la otra patada
         echo "cluster" > poscar.aux
         echo "1.0" >> poscar.aux
         cat ../programs/Matriz_$file_name >> poscar.aux
         echo $Simbolo_1 $Simbolo_2 >> poscar.aux
         echo "$N_Simbolo_1 $N_Simbolo_2" >> poscar.aux
         echo "Direct" >> poscar.aux
         cat aux >> poscar.aux
         cp ../programs/kickpp .
         cp ../programs/Matriz_$file_name .
         ./kickpp $step_width $z_vac_min Matriz_$file_name
         rm kickpp
         tail -$Nat aux2 | awk '{print $2 " "$3" "$4 }' > coords.tmp
         rm aux2
         mv coords.tmp ../programs/
         cd ../programs/
         for ((jinv=1; jinv<$(($Nat+1)); jinv++))
         do
            head -$jinv coords.tmp | tail -1 | ./inverse Matriz_$file_name >> inverted
         done
         rm coords.tmp
         mv inverted ../$file_name/preposcar
         cd ../$file_name
         rm poscar.aux
         rm Matriz_$file_name
      fi
#      rm aux

   fi                                             # Cierra el move-swap. Continua el algoritmo


   if [ $Sel -eq 1 ]
   then                                                             #Si hay selective dynamics

      for ((iauxil=0;iauxil<$Nat;iauxil++))
      do

         echo " T  T  T" >> din

      done


      tail -$Nat preposcar >> aux20
      paste aux20 din >>  POSCAR
      rm din aux20 preposcar

   else                      #De otra forma echa directamente las coordenadas de aux a POSCAR

      tail -$Nat preposcar >> POSCAR
      rm preposcar

   fi
   echo "Iteration $i of structure $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2 ($file_name)"

   ./run.sh
   contenido=$(grep "reached required accuracy" OUTCAR | wc -l )
########################################################################################## COMIENZA FASE 4

   while [[ $contenido -ne 1  ]]
   do

      echo " --> SCF failed. Starting again from randomly generated structure! "
      rm CHG CHGCAR DOSCAR EIGENVAL XDATCAR IBZKPT OSZICAR PCDAT REPORT WAVECAR *.xml CONTCAR POSCAR
      cp POSCARinitial POSCAR
cp ../programs/Matriz_$file_name .
      if [ $n -gt 3 ]
      then
        #Caso bimetálico
         if [ $randomness -eq 1 ] #aleatoriedad total corre el generador actual
         then
            ./../programs/SRandomGenerator $Simbolo_1 $N_Simbolo_1 $Simbolo_2 $N_Simbolo_2 $z_vac_min Matriz_$file_name
         else
            ./../programs/RandomGenerator $Simbolo_1 $N_Simbolo_1 $Simbolo_2 $N_Simbolo_2 $z_vac_min Matriz_$file_name
         fi
      else                                                                     #Caso monometálico
         if [ $randomness -eq 1 ] #aleatoriedad total corre el generador actual
         then
            ./../programs/SRandomGenerator $Simbolo_1 $N_Simbolo_1 $z_vac_min Matriz_$file_name
         else
            ./../programs/RandomGenerator $Simbolo_1 $N_Simbolo_1 $z_vac_min Matriz_$file_name
         fi
      fi
rm Matriz_$file_name
      tail -$(($Nat)) ClusterGenerated.xyz | awk '{print $2 "  " $3 "  " $4 }' >> aux
      rm ClusterGenerated.xyz
      if [ $direct -eq 1 ]                       #Analiza si está en formato cartesiano o directo
      then                        #Si está en direct convierte las coords del clúster a  directas
         cp aux ../programs/
         cd ../programs/
         for ((jinv=1; jinv<$(($Nat+1)); jinv++))
         do
            head -$jinv aux | tail -1 | ./inverse Matriz_$file_name >> inverted
         done
         rm aux                             #Este contiene las coordenadas Cartesianas originales
         mv inverted ../$file_name/aux                     #inverted tiene las coords directas, las manda a aux
         cd ../$file_name/
      fi


      if [ $Sel -eq 1 ]
      then                                                             #Si hay selective dynamics

         for ((iauxi=0;iauxi<$Nat;iauxi++))
         do

            echo " T T T" >> din

         done

         tail -$Nat aux >> aux2
         paste aux2 din >>  POSCAR
         rm din aux aux2

      else      #De otra forma echa directamente las coordenadas de aux a POSCAR

         tail -$Nat aux >> POSCAR
         rm aux

      fi

      ./run.sh

      contenido=$(grep "reached required " OUTCAR | wc -l )

   done #Continua con el codigo si si convergio

################################################################################################
#                                     Save  configuration                                      #
################################################################################################

   
   EnergiaAnterior=$(  head -1  CONTCAR$(($i-1)) | awk '{print $2}'   )
   Energia=$(tail -1 OSZICAR | awk '{print $5 }')    #Extrae energia del OSZICAR

   echo "$i     $Energia " >> CONTCAR$i
   N=$(wc -l CONTCAR | awk '{ print $1 }' )
   tail -$(($N-1)) CONTCAR >> CONTCAR$i  #Renombra los archivos y les agrega la energia
   mv POSCAR POSCAR$i
   rm CHG CHGCAR DOSCAR OSZICAR EIGENVAL XDATCAR IBZKPT  PCDAT REPORT WAVECAR *.xml CONTCAR

################################################################################################
#                                  Metropolis Monte-Carlo                                      #
################################################################################################

   accepted=$(python ../programs/metropolis.py $EnergiaAnterior $Energia $Temperature)

      if [ $accepted = true ]
      then   # La energía ha sido aceptada

        echo "--> Basin Hopping MC criteria: Energy accepted! "
        echo "--> Finished iteration $i"
        head -1 CONTCAR$i >> energies.txt
        tail -$i energies.txt |  sort -nk2 > sorted.txt
        i=$(($i+1)) #Convergencia lograda, aumenta en 1 el contador

      else

        echo "--> Basin Hopping MC criteria: Energy rejected!"
        mv CONTCAR$i rejected/CONTCARrejected$i
        mv POSCAR$i rejected/POSCARrejected$i
        m=$(($m+1))

      fi


   done

cd ../programs
rm Matriz_$file_name
cd ..

cd $file_name
echo "About this project: $file_name" >> Summary.txt
echo "
input.bh file used
============================input.bh=============================
" >>Summary.txt
cat ../input.bh >> Summary.txt
echo "
=================================================================
" >>Summary.txt
echo "
INCAR file used
=============================INCAR===============================
" >> Summary.txt

cat  ../input/INCAR >> Summary.txt
echo "
=================================================================
KPOINTS file used
=============================KPOINTS==============================
">> Summary.txt

cat ../input/KPOINTS >> Summary.txt

echo "
==================================================================" >> Summary.txt
