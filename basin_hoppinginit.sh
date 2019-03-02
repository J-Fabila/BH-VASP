#! /bin/bash

################################################################################################
#                                    Get data from input.bh                                    #
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
NPOSCAR=$(cat input/POSCAR | grep . | wc -l ) #Numero de lineas del poscar sin cluster


##############################################################################################
#                                         BEGIN ALGORITHM                                    #
##############################################################################################

if [ $NPOSCAR -gt 5 ]                        #Determina si el sistema es gas phase o soportado
then                   #Si es soportado agrega los simbolos atomicos del cluster y sus numeros
   cd input
   cat POSCAR | grep . >> aux
   rm POSCAR
   cat aux >> POSCAR
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
else          #Si es gas phase prepara el POSCAR, agrega los simbolos atomicos y el "Cartesian"
   cd input
   cat POSCAR | grep . >> aux
   rm POSCAR
   cat aux >> POSCAR
   rm aux

   echo " $Simbolo_1 $Simbolo_2" >> POSCAR
   echo "$N_Simbolo_1 $N_Simbolo_2 " >> POSCAR
   echo "Cartesian">>POSCAR

fi


################################################################################################
#                                   Generates POTCAR file                                      #
################################################################################################



head -6 POSCAR | tail -1 >> Species                #Toma la linea 6 del POSCAR, que contiene los
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


################################ GENERATES RUN.sh FILE #######################################

echo "mpirun -np $Ncore vasp_gam > output.out" >> run.sh
chmod +x run.sh

########## Takes the matrix from POSCAR file, this is necessary for calculations #############

head -5 POSCAR | tail -3 | awk '{print $1 " " $2 "  " $3}' >> ../programs/Matriz

##############################################################################################

cd ..

cp -r input $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2  ## OJO ACA: REVISAR COMO COPIA #Copia los archivos de input al directorio de trabajo
#python programs/RandomGenerator.py $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2/POSCAR   #Genera una estructura aleatoria y concatena a POSCAR

cd $Simbolo_1$N_Simbolo_1$Simbolo_2$N_Simbolo_2         #Nos mueve a ese directorio de trabajo

echo "  " >>aux      #Genera el archivo aux, ahi se pondrán las coordenadas de RandomGenerator


if [ $Npath -gt 8 ]                           #Determina si hay archivo de inicializacion o no
then                       #Si existe toma las coordenadas del mismo y las concatena al POSCAR

   tail -$Nat $path >>  POSCAR

else                                             #De otra forma  invoca al generador aleatorio
      if [ $n -gt 3 ] #Determina y corre de acuerdo con si es bimetálico  o monometálico #El 3 s arbitrario, solo determina si el segundo campo es vacio o no
      then                                                                    #Caso bimetalico
         python ../programs/RandomGenerator.py aux $Nt1,$Nt2 $XRange $YRange $ZRange $ZVacuum
      else                                                                  #Caso monometalico
         python ../programs/RandomGenerator.py aux $Nt1 $XRange $YRange $ZRange $ZVacuum
      fi
fi

                              #Si hay selective dynamics activado entonces le coloca "T" a las
                                         # coordenadas del clúster, de otra forma los deja así

if [ $Sel -eq 1 ] 
then                                                                #Si hay selective dynamics

   for ((i=0;i<$Nat;i++))
   do

      echo " T  T  T">>din

   done


   tail -$Nat aux >> aux2
   paste aux2 din >>  POSCAR

   rm din aux aux2

else                          #De otra forma echa directamente las coordenadas de aux a POSCAR 

   tail -$Nat aux >> POSCAR
fi   

mkdir rejected                                               #Crea el directorio de rechazados
./run.sh                                                    #Crea el subproceso que corre VASP
contenido=$(grep "reached required " OUTCAR | wc -l )        # 0 si no converge, 1 si converge


#ACA INICIA LA  FASE 2 DE PRUEBA
while [[ $contenido -ne 1  ]]       #Mientras no converja la configuracion reinicia el calculo
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

   if [ $Sel -eq 1 ] 
   then                                                             #Si hay selective dynamics

      for ((i=0;i<$Nat;i++))
      do

         echo " T  T  T">>din

      done


      tail -$Nat aux >> aux2
      paste aux2 din >>  POSCAR

      rm din aux aux2

   else                      #De otra forma echa directamente las coordenadas de aux a POSCAR 

      tail -$Nat aux >> POSCAR
   fi

   ./run.sh
   contenido=$(grep "reached required " OUTCAR | wc -l )


done                                                   #Continua con el codigo si si convergio


Energia=$(tail -1 OSZICAR | awk '{print $5 }')                 #Extrae  la energia del OSZICAR
echo "1         $Energia " >> CONTCAR1      #Escribe la energia y num de iteracion en CONTCAR1
N=$(wc -l CONTCAR | awk '{ print $1 }' )         #Cuenta el numero de lineas que tiene CONTCAR
tail -$(($N-1)) CONTCAR >> CONTCAR1       #Mueve la informacion a CONTCAR1 con el nuevo titulo
mv POSCAR POSCAR1                                                      #Mueve POSCAR a POSCAR1
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

#COMIENZA FASE 3 DE PRUEBA

i=2

while [ $i -lt $(($iteraciones+1)) ]
do

   let resto=$i%10                          # Calcula el resto de i  para hacer los pasos swap

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
   then                                                             #Si hay selective dynamics

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

   while [[ $contenido -ne 1  ]]
   do
      echo " --> SCF failed. Starting again from randomly generated structure! "

      if [ $Sel -eq 1 ]              # Este if extrae las coordenadas de la iteracion anterior
      then                                       # dependiendo del caso: soportado o gas phase

         head -$(($NPOSCAR+$Nat)) CONTCAR$(($i+1)) >> aux2
         echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR
         head -$NPOSCAR CONTCAR1 | tail -$(($NPOSCAR-1)) >> POSCAR
         tail -$Nat aux2 >> aux              #Este contiene las coordenadas que leerá despues
         rm aux2

      else

         head -$(($NPOSCAR+3+$Nat)) CONTCAR$(($i-1)) >> aux2
         echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR
         head -8 CONTCAR1 | tail -7 >> POSCAR
         tail -$Nat aux2 >> aux              #Este contiene las coordenadas que leerá despues

      fi

      rm CHG CHGCAR DOSCAR EIGENVAL XDATCAR IBZKPT OSZICAR PCDAT REPORT WAVECAR *.xml CONTCAR POSCAR

      echo "BH_DFT_VASP: POSCAR of iteration $i" >> POSCAR

      cat aux >> POSCAR
      rm aux

      echo "   " >> auxtoinvert            #Crea el archivo donde se guardarán las coordenadas

      if [ $n -gt 3 ] #Determina y corre de acuerdo con si es bimetálico  o monometálico
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

