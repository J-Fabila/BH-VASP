#include<atomicpp.h>

int main(int argc, char *argv[])
{
   Cluster np;
   np.read_VASP(argv[1]);
   np.print_xyz(argv[2]);
   return 0;
}
