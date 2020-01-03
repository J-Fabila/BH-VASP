#include"atomicpp.h"

int main(int argc, char *argv[])
{
   double xm, ym, zm;
   double dx,dy,dz;
   double matrix[3][3];
   Cluster clus;
   if(argc==5)
   {
      clus.rand_generator(argv[1],atoi(argv[2]));
      zm=stod(argv[3]);
      ifstream matriz(argv[4]);
      for(int i=0;i<3;i++)
      {
         matriz>>matrix[i][0]>>matrix[i][1]>>matrix[i][2];
      }
   }
   else
   {
      if(argc==7)
      {
         clus.rand_generator(argv[1],atoi(argv[2]),argv[3],atoi(argv[4]));
         zm=stod(argv[5]);
         ifstream matriz(argv[6]);
         for(int i=0;i<3;i++)
         {
            matriz>>matrix[i][0]>>matrix[i][1]>>matrix[i][2];
         }
      }
   }
   dx=(matrix[0][0]+matrix[1][0]+matrix[2][0])/3.0;
   dy=(matrix[0][1]+matrix[1][1]+matrix[2][1])/3.0;
   double zmin=clus.z_min();

   dx=dx+random_number(-2.0,2.0);
   dy=dy+random_number(-2.0,2.0);
   dz=-zmin+zm;

   clus.move(dx,dy,dz);
   clus.print_xyz("ClusterGenerated.xyz");
   return 0;
}

