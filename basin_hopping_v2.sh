################################################################################################
#    Basin Hopping Algorithm coupled with VASP
#
#    Author:
#    Jorge Fabila
#    https://github.com/J-Fabila
#
#    Oliver Paz-Borbón
#
################################################################################################



################################################################################################
#                                    Gets data from input.bh                                   #
################################################################################################

Simbolo_1=$(grep "cluster_ntyp" input.bh | cut -d "[" -f 2 | cut -d ":" -f 1)
Simbolo_2=$(grep "cluster_ntyp" input.bh | cut -d "[" -f 3 | cut -d ":" -f 1) 2>/dev/null
N_Simbolo_1=$(grep "cluster_ntyp" input.bh | cut -d "[" -f 2 | cut -d ":" -f 2 | cut -d "]" -f 1)
N_Simbolo_2=$(grep "cluster_ntyp" input.bh | cut -d "[" -f 3 | cut -d ":" -f 2 | cut -d "]" -f 1) 2>/dev/null
Nt1=$(grep "cluster_ntyp"  input.bh | cut -d " " -f 3 | cut -d "," -f 1 )
Nt2=$(grep "cluster_ntyp"  input.bh | cut -d " " -f 4) 2>/dev/null
#Alternativa, a ver si funcioa
#Nt1=$(grep "cluster_ntyp"  input.bh | cut -d "=" -f 2 | cut -d "," -f 1 )
#Nt2=$(grep "cluster_ntyp"  input.bh | cut -d "," -f 2)

n=$(echo $Nt2  | wc -c  )  #Este es un criterio para determinar si es bimetálico o no
XRange=$(grep "x_range"  input.bh | cut -d " " -f 3)
YRange=$(grep "y_range"  input.bh | cut -d " " -f 3)
ZRange=$(grep "z_range"  input.bh | cut -d " " -f 3)
ZVacuum=$(grep "z_vacuum"  input.bh | cut -d " " -f 3)

xmin=$(grep "x_range" input.bh | cut -d "[" -f 2 | cut -d ":" -f 1 )
xmax=$(grep "x_range" input.bh | cut -d "[" -f 2 | cut -d ":" -f 2 | cut -d "]" -f 1)

ymin=$(grep "y_range" input.bh | cut -d "[" -f 2 | cut -d ":" -f 1 )
ymax=$(grep "y_range" input.bh | cut -d "[" -f 2 | cut -d ":" -f 2 | cut -d "]" -f 1)

zmin=$(grep "z_range" input.bh | cut -d "[" -f 2 | cut -d ":" -f 1 )
zmax=$(grep "z_range" input.bh | cut -d "[" -f 2 | cut -d ":" -f 2 | cut -d "]" -f 1)

z_vac_min=$(grep "z_vacuum" input.bh | cut -d "[" -f 2 | cut -d ":" -f 1 )
z_vac_max=$(grep "z_vacuum" input.bh | cut -d "[" -f 2 | cut -d ":" -f 2 | cut -d "]" -f 1)


#############################################################################
randomness= $(grep "randomness" input.bh | awk '{print $3}')
#############################################################################

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
Ncore=$(grep "Ncore" input.bh | awk '{print $3}')
iteraciones=$(grep "iterations" input.bh | awk '{ print $3 }' )
path=$(grep "initialization_file" input.bh | awk '{ print $3 }' )
Npath=$(echo $path | wc -c )
m=1 #contador de coords rechazadas
Sel=$(grep "Selective"  input/POSCAR | wc -l )   #Determina si hay un selective dynamics
NPOSCAR=$(cat input/POSCAR | grep . | wc -l ) #Numero de lineas del poscar sin cluster


swap_step=10      # SWAP STEP

##############################################################################################
#                                         BEGIN ALGORITHM                                    #
##############################################################################################

cd programs
rm Matriz 2> /dev/null
cd ..

if [ $NPOSCAR -gt 8 ]                        #Determina si el sistema es gas phase o soportado
then                   #Si es soportado agrega los simbolos atomicos del cluster y sus numeros
   cd input
   cat POSCAR | grep . >> aux
   rm POSCAR
   cat aux >> POSCAR
   rm aux
                                                                #-----------------------------------#
   Species=$(head -6 POSCAR | tail -1)                          #                                   #
   Numbers=$(head -7 POSCAR | tail -1)                          #                                   #
                                                                #            Estas lineas           #
   head -5 POSCAR >> POSCARinitial                              #       concatenan los Símbolos     #
   echo "$Species $Simbolo_1 $Simbolo_2" >> POSCARinitial       #            y números del          #
   echo "$Numbers $N_Simbolo_1 $N_Simbolo_2 " >> POSCARinitial  #         clúster al POSCAR         #
   tail -$(($NPOSCAR-7)) aux >> POSCARinitial                   #                                   #
   rm aux                                                       #                                   #
                                                                #-----------------------------------#
else          #Si es gas phase prepara el POSCAR, agrega los simbolos atomicos y el "Cartesian"
   cd input
   cat POSCAR | grep . >> aux
   rm POSCAR
   cat aux >> POSCAR
   rm aux

   head -5 POSCAR >> POSCARinitial
   echo " $Simbolo_1 $Simbolo_2" >> POSCARinitial
   echo "$N_Simbolo_1 $N_Simbolo_2 " >> POSCARinitial
   echo "Cartesian">>POSCARinitial

fi

direct=$(grep "irect" POSCAR | wc -l )

################################################################################################
#                                   Generates POTCAR file                                      #
################################################################################################


head -6 POSCARinitial | tail -1 >> Species                #Toma la linea 6 del POSCAR, que contiene los
                                        # simbolos atomicos y los guarda en el archivo "Species"
tr -s '[:blank:]' '\n' < Species >> aux  # Convierte la linea 6 en una columna, lo guarda en aux
grep . aux >> Symbols         #Quita espacios en blanco (lineas vacias) y lo guarda en "Symbols"
rm aux                                                          # Ya no se necesita este archivo
Ntyp=$(cat Symbols | wc -l ) #El numero de simbolos (Ntyp) es igual al numero de lineas ocupadas
rm Species
#                Cada columna tiene un simbolo atomico, lee cada una y lo guarda como variable S

for((i=1;((i<$Ntyp+1));i++))
do
   S=$(head -$i Symbols |tail -1 )
   cat $Pseudo_Dir/$S$pseudotype/POTCAR>>POTCAR          #Concatena desde el PseudoDir al POTCAR
done

rm Symbols


################################## GENERATES RUN.sh FILE #######################################

echo "mpirun -np $Ncore vasp_gam > output.out" >> run.sh
chmod +x run.sh

########## Takes the matrix from POSCAR file, this is necessary for calculations #############

head -5 POSCAR | tail -3 | awk '{print $1 " " $2 "  " $3}' >> ../programs/Matriz

##############################################################################################

cd ..
################################# CREATES Work Directory #####################################

if [ -d $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2 ]
then
   mv $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2 $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2_other
fi

cp -r input $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2


###############################################################################################
rm input/POSCARinitial									      #
rm input/POTCAR
rm input/run.sh
cd $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2         #Nos mueve a ese directorio de trabajo#
rm POSCAR										      #
cp POSCARinitial POSCAR									      #
#####################################################################Auxiliar para no dejar el poscar modificado

if [ $Npath -gt 8 ]                           #Determina si hay archivo de inicializacion o no
then                           #Si existe toma las coordenadas del mismo y las copia al POSCAR

   cp $path POSCAR

else                                             #De otra forma  invoca al generador aleatorio

   if [ $n -gt 3 ]          #Determina y corre de acuerdo con si es bimetálico  o monometálico
   then                                                                       #Caso bimetálico
      if [ $randomness -eq 1 ] #aleatoriedad total corre el generador actual
      then
         ./../programs/SRandomGenerator $Simbolo_1 $N_Simbolo_1 $Simbolo_2 $N_Simbolo_2 $xmin $xmax $ymin $ymax $zmin $zmax $z_vac_min $z_vac_max
      else
         ./../programs/RandomGenerator $Simbolo_1 $N_Simbolo_1 $Simbolo_2 $N_Simbolo_2 $xmin $xmax $ymin $ymax $zmin $zmax $z_vac_min $z_vac_max
      fi
   else                                                                     #Caso monometálico
      if [ $randomness -eq 1 ] #aleatoriedad total corre el generador actual
      then
         ./../programs/SRandomGenerator $Simbolo_1 $N_Simbolo_1 $xmin $xmax $ymin $ymax $zmin $zmax $z_vac_min $z_vac_max
      else
         ./../programs/RandomGenerator $Simbolo_1 $N_Simbolo_1 $xmin $xmax $ymin $ymax $zmin $zmax $z_vac_min $z_vac_max
      fi
   fi

   tail -$(($Nat+1)) ClusterGenerated.xyz | awk '{print $2 "  " $3 "  " $4 }' >> aux
   rm ClusterGenerated.xyz

   if [ $direct -eq 1 ]                       #Analiza si está en formato cartesiano o directo
   then                        #Si está en direct convierte las coords del clúster a  directas

      for ((jinv=1; jinv<$(($Nat+1)); jinv++))
      do

         cd ../programs ; head -$jinv aux | tail -1 | ./inverse >> ../$Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2/inverted
           cd ../$Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2

      done
      rm aux                             #Este contiene las coordenadas Cartesianas originales
      mv inverted aux                     #inverted tiene las coords directas, las manda a aux

   fi

fi

                              #Si hay selective dynamics activado entonces le coloca "T" a las
                                         # coordenadas del clúster, de otra forma los deja así

if [ $Sel -eq 1 ]
then                                                                #Si hay selective dynamics

   for ((iaux=0;iaux<$Nat;iaux++))
   do

      echo " T  T  T">>din

   done


   tail -$Nat aux >> aux2
   paste aux2 din >>  POSCAR

   rm din aux aux2

else                          #De otra forma echa directamente las coordenadas de aux a POSCAR

   tail -$Nat aux >> POSCAR
   rm aux
fi

mkdir rejected                                               #Crea el directorio de rechazados
cat POSCAR
./run.sh                                                    #Crea el subproceso que corre VASP
contenido=$(grep "reached required " OUTCAR | wc -l )        # 0 si no converge, 1 si converge
echo "Configuracion CONTCAR!, terminando fase 1"
cat CONTCAR  #OJO BORRAR DESPUES

#############################################################ACA INICIA LA  FASE 2 DE PRUEBA



while [[ $contenido -ne 1  ]]       #Mientras no converja la configuracion reinicia el calculo
do

   echo " --> SCF failed. Starting again from randomly generated structure! "
   rm CHG CHGCAR DOSCAR EIGENVAL XDATCAR IBZKPT OSZICAR PCDAT REPORT WAVECAR *.xml CONTCAR POSCAR
   cp POSCARinitial POSCAR
   echo "  " >>aux

      if [ $n -gt 3 ]          #Determina y corre de acuerdo con si es bimetálico  o monometálico
      then                                                                       #Caso bimetálico
         if [ $randomness -eq 1 ] #aleatoriedad total corre el generador actual
         then
            ./../programs/SRandomGenerator $Simbolo_1 $N_Simbolo_1 $Simbolo_2 $N_Simbolo_2 $xmin $xmax $ymin $ymax $zmin $zmax $z_vac_min $z_vac_max
         else
            ./../programs/RandomGenerator $Simbolo_1 $N_Simbolo_1 $Simbolo_2 $N_Simbolo_2 $xmin $xmax $ymin $ymax $zmin $zmax $z_vac_min $z_vac_max
         fi
      else                                                                     #Caso monometálico
         if [ $randomness -eq 1 ] #aleatoriedad total corre el generador actual
         then
            ./../programs/SRandomGenerator $Simbolo_1 $N_Simbolo_1 $xmin $xmax $ymin $ymax $zmin $zmax $z_vac_min $z_vac_max
         else
            ./../programs/RandomGenerator $Simbolo_1 $N_Simbolo_1 $xmin $xmax $ymin $ymax $zmin $zmax $z_vac_min $z_vac_max
         fi
      fi

   tail -$(($Nat+1)) ClusterGenerated.xyz | awk '{print $2 "  " $3 "  " $4 }' >> aux
   rm ClusterGenerated.xyz

   if [ $direct -eq 1 ]                       #Analiza si está en formato cartesiano o directo
   then                        #Si está en direct convierte las coords del clúster a  directas

      for ((jinv=1; jinv<$(($Nat+1)); jinv++))
      do

         cd ../programs ; head -$jinv aux | tail -1 | ./inverse >> ../$Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2/inverted
           cd ../$Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2

      done
      rm aux                             #Este contiene las coordenadas Cartesianas originales
      mv inverted aux                     #inverted tiene las coords directas, las manda a aux

   fi



   if [ $Sel -eq 1 ]
   then                                                             #Si hay selective dynamics

      for ((iauxi=0;iauxi<$Nat;iauxi++))
      do

         echo " T  T  T">>din

      done


      tail -$Nat aux >> aux2
      paste aux2 din >>  POSCAR

      rm din aux aux2

   else                      #De otra forma echa directamente las coordenadas de aux a POSCAR

      tail -$Nat aux >> POSCAR
      rm aux

   fi

echo "POSCAR CONFIGURACION"  #OJO ACA:BORRA ESTS LINEAS CUAND ACABES Con las pruebas
cat POSCAR
   ./run.sh
echo "CONTCAR RELAJADO"
cat CONTCAR

   contenido=$(grep "reached required " OUTCAR | wc -l )


done                                                   #Continua con el codigo si si convergio


Energia=$(tail -1 OSZICAR | awk '{print $5 }')                 #Extrae  la energia del OSZICAR
echo "1     $Energia " >> CONTCAR1      #Escribe la energia y num de iteracion en CONTCAR1
N=$(wc -l CONTCAR | awk '{ print $1 }' )         #Cuenta el numero de lineas que tiene CONTCAR
tail -$(($N-1)) CONTCAR >> CONTCAR1       #Mueve la informacion a CONTCAR1 con el nuevo titulo
mv POSCAR POSCAR1                                                      #Mueve POSCAR a POSCAR1
cat CONTCAR1   #OJOACA: AUXILIAR BORRARDESPUES
rm CHG CHGCAR DOSCAR EIGENVAL XDATCAR IBZKPT OSZICAR PCDAT REPORT WAVECAR *.xml CONTCAR

echo "Step  Energy" >> energies.txt
echo "1     $Energia" >> energies.txt

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

##########################################################################COMIENZA FASE 3 DE PRUEBA

i=2

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

         cd ../programs ; echo "$dx $dy $dz" | ./inverse > ../$Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2/kick
echo "=====================kick========================="
cat ../$Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2/kick
echo "=================================================="
         cd ../$Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2


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
         for ((j=1; j<$(($N_Simbolo_1+1)); j++))
         do
            echo $Simbolo_1 >> simbolos.aux
         done
         for ((j=1; j<$(($N_Simbolo_2+1)); j++))
         do
            echo $Simbolo_2 >> simbolos.aux
         done
         paste simbolos.aux aux >> xyz
         echo $Nat >> cluster.xyz
         echo " " >> cluster.xyz
         cat xyz >> cluster.xyz
         #kickpp.cpp lee cluster.xyz.Hace la patada
         cp ../programs/kickpp .
         ./kickpp # < cluster.xyz > aux2
         rm kickpp
         tail -$Nat aux2 >> preposcar                 #Este contiene las coordenadas que leerá despues
         rm cluster.xyz
         rm xyz
         rm simbolos.aux
      fi
      rm aux

   fi                                             # Cierra el move-swap. Continua el algoritmo


   if [ $Sel -eq 1 ]
   then                                                             #Si hay selective dynamics

      for ((iauxil=0;iauxil<$Nat;iauxil++))
      do

         echo " T  T  T" >> din

      done


      tail -$Nat preposcar >> aux20
      paste aux20 din >>  POSCAR
      echo "aux20"
      cat aux20
      echo "din"
      cat din
      rm din aux20 preposcar

   else                      #De otra forma echa directamente las coordenadas de aux a POSCAR

      tail -$Nat preposcar >> POSCAR
      rm preposcar

   fi
  cat POSCAR
   echo "Iteration $i of structure $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2"

   ./run.sh
   echo "CONTCAR terminando fase 3"
   cat CONTCAR  #OJO: AUXILIAR BORRAR DESPUES
   contenido=$(grep "reached required accuracy" OUTCAR | wc -l )
########################################################################################## COMIENZA FASE 4

   while [[ $contenido -ne 1  ]]
   do

      echo " --> SCF failed. Starting again from randomly generated structure! "
      rm CHG CHGCAR DOSCAR EIGENVAL XDATCAR IBZKPT OSZICAR PCDAT REPORT WAVECAR *.xml CONTCAR POSCAR
      cp POSCARinitial POSCAR

      if [ $n -gt 3 ]
      then
        #Caso bimetálico
         if [ $randomness -eq 1 ] #aleatoriedad total corre el generador actual
         then
            ./../programs/SRandomGenerator $Simbolo_1 $N_Simbolo_1 $Simbolo_2 $N_Simbolo_2 $xmin $xmax $ymin $ymax $zmin $zmax $z_vac_min $z_vac_max
         else
            ./../programs/RandomGenerator $Simbolo_1 $N_Simbolo_1 $Simbolo_2 $N_Simbolo_2 $xmin $xmax $ymin $ymax $zmin $zmax $z_vac_min $z_vac_max
         fi
      else                                                                     #Caso monometálico
         if [ $randomness -eq 1 ] #aleatoriedad total corre el generador actual
         then
            ./../programs/SRandomGenerator $Simbolo_1 $N_Simbolo_1 $xmin $xmax $ymin $ymax $zmin $zmax $z_vac_min $z_vac_max
         else
            ./../programs/RandomGenerator $Simbolo_1 $N_Simbolo_1 $xmin $xmax $ymin $ymax $zmin $zmax $z_vac_min $z_vac_max
         fi
      fi
      tail -$(($Nat+1)) ClusterGenerated.xyz | awk '{print $2 "  " $3 "  " $4 }' >> aux
      rm ClusterGenerated.xyz

      if [ $direct -eq 1 ]                       #Analiza si está en formato cartesiano o directo
      then                        #Si está en direct convierte las coords del clúster a  directas

         for ((jinv=1; jinv<$(($Nat+1)); jinv++))
         do

            cd ../programs ; head -$jinv aux | tail -1 | ./inverse >> ../$Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2/inverted
            cd ../$Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2

         done
         rm aux                             #Este contiene las coordenadas Cartesianas originales
         mv inverted aux                     #inverted tiene las coords directas, las manda a aux

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
echo "POSCAR de random generator, fase 4"
cat POSCAR  #OJO: AXULIAR NOMAS; BORRAR LUego

      ./run.sh
echo "CONTCAR de random generator, fase 4"
cat CONTCAR

      contenido=$(grep "reached required " OUTCAR | wc -l )

   done #Continua con el codigo si si convergio

################################################################################################
#                                     Save  configuration                                      #
################################################################################################

   EnergiaAnterior=$(echo $Energia)
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
rm Matriz
cd ..

cd $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2
echo "About this project:" >> Summary.txt
echo "
INCAR file used
==============================INCAR===============================
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
