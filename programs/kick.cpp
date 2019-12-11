#include<atomicpp.h>

int main()
{
   Cluster Np;
   Np.read_VASP("CONTCAR");
   Np.centroid();
   Np.kick_lennard(1.5);
   Np.print_xyz("configuration.xyz");
   //system("cat Matriz >> POSCAR");
   //system(" ...etc");
   //system("cat configuration.xyz >> POSCAR");
   //system("rm configuration.xyz");
   return 0;
}
