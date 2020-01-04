#Download files
echo "

"
echo -e "\e[34m Downloading files ... \e[0m"
echo "

"
sleep 1.5
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/input.bh
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/basin_hopping.sh
chmod +x basin_hopping.sh
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/continue.sh
chmod +x continue.sh
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/queue.sh
chmod +x queue.sh
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/Instruction_Manual.pdf
#Download programs
echo "

"
echo -e "\e[34m Creating directories ... \e[0m"
echo "

"
sleep 1.5
mkdir programs
cd programs
#Download files
echo "

"
echo -e "\e[34m Downloading source code ... \e[0m"
echo "

"
sleep 1.5
wget https://raw.githubusercontent.com/J-Fabila/Atomic.h/master/atomicpp.h
#Tal vez sea mejor descargar diretamente del repositorio de atomicc.pp.h
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/programs/read_crystal.cpp
g++ -o read_crystal read_crystal.cpp
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/programs/read_matriz.cpp
g++ -o read_matriz read_matriz.cpp
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/programs/kickpp.cpp
g++ -o kickpp kickpp.cpp
##Tal vez sea mejor pasar metropolis.py a cpp
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/programs/metropolis.py
## Idem
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/programs/move.py
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/programs/RandomGenerator.cpp
g++ -o RandomGenerator RandomGenerator.cpp
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/programs/SRandomGenerator.cpp
g++ -o SRandomGenerator SRandomGenerator.cpp
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/programs/inverse.c
gcc -o inverse inverse.c -lm
#A lo mejor sea mejor descargarlos directo del xyz_to_vasp
echo "

"
echo -e "\e[34m Compiled ... \e[0m"
echo "

"
sleep 1.7
mkdir auxiliar_programs
cd auxiliar_programs

wget https://raw.githubusercontent.com/J-Fabila/VASP_to_XYZ/master/XYZtoPOSCAR.sh
wget https://raw.githubusercontent.com/J-Fabila/VASP_to_XYZ/master/POSCARtoXYZ.c
mv XYZtoPOSCAR.sh to_vasp
chmod +x to_vasp
gcc -o to_xyz POSCARtoXYZ.c
#Agrega estos programas como comandos cualquiera
echo "

"
echo -e "\e[34m Installing auxiliar programs ... \e[0m"
echo "

"
sleep 1.7


dir=$PWD
cd
echo 'PATH=$PATH:'$dir >>.bashrc
source .bashrc
cd $dir
#Download input files
cd ../..
mkdir input
cd input
echo "

"
echo -e "\e[34m Downloading input files ... \e[0m"
echo "

"
sleep 1.3
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/input/INCAR
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/input/KPOINTS
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/input/POSCAR
cd ..

rm $0

echo -e "


\e[34m                     Done!


\e[0m"
sleep 1.8
clear
