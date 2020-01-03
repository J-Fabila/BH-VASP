#include<iostream>
#include <fstream>
double matrix[3][3];
using namespace std;
int main(int argc, char *argv[])
{
ifstream matriz(argv[1]);
for(int i=0;i<3;i++)
{
matriz>>matrix[i][0]>>matrix[i][1]>>matrix[i][2];
}
cout<<"0.0 "<<matrix[0][0]<<" 0.0 "<<matrix[1][1]<<" 0.0 "<<matrix[2][2]<<" 0.0 "<<matrix[2][2]<<endl;
return 0;
}
