#include"atomicpp.h" 

int main(int argc, char *argv[])
{
double _x, _y, _z;
Cluster clus;
if(argc==11)
{
  clus.srand_generator(argv[1],atoi(argv[2]));
  bool acep=clus.fit_in(atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[9]),atoi(argv[10]));
  while (acep==false)
  {
     clus.~Cluster();
     clus.srand_generator(argv[1],atoi(argv[2]));
     acep=clus.fit_in(atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[9]),atoi(argv[10]));
   }
  _x= atoi(argv[3]) + ((double)rand())/((double)RAND_MAX )* (atoi(argv[4])-atoi(argv[3]));
  _y= atoi(argv[5]) + ((double)rand())/((double)RAND_MAX )* (atoi(argv[6])-atoi(argv[5]));
  _z= atoi(argv[9]);

}
else
{
if(argc==13)
{

clus.srand_generator(argv[1],atoi(argv[2]),argv[3],atoi(argv[4]));
bool acep=clus.fit_in(atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),atoi(argv[11]),atoi(argv[12]));
while (acep==false)
{
   clus.~Cluster();
   clus.srand_generator(argv[1],atoi(argv[2]),argv[3],atoi(argv[4]));
   acep=clus.fit_in(atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),atoi(argv[11]),atoi(argv[12]));
 }
_x= atoi(argv[5]) + ((double)rand())/((double)RAND_MAX )* (atoi(argv[6])-atoi(argv[5]));
_y= atoi(argv[7]) + ((double)rand())/((double)RAND_MAX )* (atoi(argv[8])-atoi(argv[7]));
_z= atoi(argv[11]);
}
}
double zmin=clus.z_min();
clus.move( _x , _y , -zmin + _z );

  clus.print_xyz("ClusterGenerated.xyz");
  float rnd=random_number(10,20);
  return 0;
}
