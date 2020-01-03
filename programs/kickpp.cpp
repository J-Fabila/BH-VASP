#include"atomicpp.h"
int main(int argc, char *argv[])
{
   double matrix[3][3];
   ifstream matriz(argv[3]);
   for(int i=0;i<3;i++)
   {
      matriz>>matrix[i][0]>>matrix[i][1]>>matrix[i][2];
   }

   Cluster Np;
   double dx,dy,dz,zmin,zm;
   Np.read_VASP("poscar.aux");
   Np.centroid();
   Np.kick_lennard(stof(argv[1]));

   dx=(matrix[0][0]+matrix[1][0]+matrix[2][0])/3.0;
   dy=(matrix[0][1]+matrix[1][1]+matrix[2][1])/3.0;
   zmin=Np.z_min();
   zm=stod(argv[2]);

   dx=dx+random_number(-2.0,2.0);
   dy=dx+random_number(-2.0,2.0);
   dz=-zmin+zm;

   Np.move(dx,dy,dz);
   Np.print_xyz("aux2");
   return 0;
}
