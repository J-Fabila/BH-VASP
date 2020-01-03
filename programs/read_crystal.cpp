#include"atomicpp.h"
int main()
{
   Crystal cris;
   double Dz;
   double xmin,xmax,ymin,ymax,zmin,zmax,zvacmin,zvacmax;
   cris.read_VASP("POSCAR");
//   Dz=cris.z_min();
//Dz=(-1.0)*Dz;
  // cris.move(0.0,0.0,Dz);
   xmin=cris.x_min();
   xmax=cris.x_max();
   ymin=cris.y_min();
   ymax=cris.y_max();
   zmin=0.0;
   zmax=cris.lattice[2][2];
   zvacmin=cris.z_max();
   zvacmax=zmax+cris.z_min();
   cout<<xmin<<" "<<xmax<<" "<<ymin<<" "<<ymax<<" "<<zmin<<" "<<zmax<<" "<<zvacmin<<" "<<zvacmax<<endl;
   return 0;
}
