#include<stdio.h>
#include<stdlib.h>
#include<math.h>

float k[3];
float kt[3];
float M[3][3];
float Mt[3][3];
float s;


int i;
int j;
float det;
int c, m, n;
float M[3][3];
float aux[3][3];
float inverse[3][3];

int main(int argc, char *argv[])
{
scanf("%f  %f  %f\n", &k[0], &k[1], &k[2]);

FILE *f = fopen(argv[1], "r");

for(i=0;i<3;i++) //Lee los elementos de matriz
{
   fscanf(f,"%f %f %f \n", &Mt[i][0], &Mt[i][1], &Mt[i][2]);
//   printf(" %f %f %f \n", Mt[i][0], Mt[i][1], Mt[i][2]);
}

for(i=0;i<3;i++)
{
for(j=0;j<3;j++)
{
M[i][j]=Mt[j][i];
}
}


det=(M[0][0]*((M[1][1]*M[2][2])-(M[1][2]*M[2][1])))-(M[0][1]*((M[1][0]*M[2][2])-(M[1][2]*M[2][0])))+(M[0][2]*((M[1][0]*M[2][1])-(M[1][1]*M[2][0])));

//printf("%f\n",det);

aux[0][0]=((M[1][1]*M[2][2])-(M[2][1]*M[1][2]));
aux[0][1]=((M[1][0]*M[2][2])-(M[2][0]*M[1][2]));
aux[0][2]=((M[1][0]*M[2][1])-(M[2][0]*M[1][1]));
aux[1][0]=((M[0][1]*M[2][2])-(M[2][1]*M[0][2]));
aux[1][1]=((M[0][0]*M[2][2])-(M[2][0]*M[0][2]));
aux[1][2]=((M[0][0]*M[2][1])-(M[2][0]*M[0][1]));
aux[2][0]=((M[0][1]*M[1][2])-(M[1][1]*M[0][2]));
aux[2][1]=((M[0][0]*M[1][2])-(M[1][0]*M[0][2]));
aux[2][2]=((M[0][0]*M[1][1])-(M[1][0]*M[0][1]));


for(i=0;i<3;i++) //Transpone, multiplica por el determinante y saca cofactores/levi civita)
{
for(j=0;j<3;j++)
{
inverse[i][j]=(1/det)*pow(-1,i+j)*aux[j][i];
//printf("%f\t",inverse[i][j]);
}
//printf("\n");
}
for(i=0;i<3;i++) //Aplica el operador al vector kick
{
  s=0;
  for(j=0;j<3;j++)
  {
  s=k[j]*inverse[i][j]+s;
  }
  kt[i]=s;
}

printf("%f %f %f\n", kt[0], kt[1], kt[2]);
fclose(f);


return 0;
}


