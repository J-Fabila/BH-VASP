//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_ RandomGenerator.cpp _/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

#include<iostream>
#include<stdlib.h>
#include<string.h>
#include <cstdlib>
#include <fstream>
#include <stdio.h>
#include  <string>
#include <sstream>
#include   <cmath>
#include   <ctime>

#include <string>
#include <utility>
#include <map>
#include <iomanip>

using namespace std;

/********* * * * * * *  *  *   *   *   *  *  *  * * * ***************/

int Nat;
int i,j,k;

// Patch for older compilers 
namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

/********************************************************************/
/************************* Atom Definition **************************/
/********************************************************************/

class Atom
{
   public:
     string Symbol;
     double x[3];
     double R;
     Atom();
     void read_Atom(string, float, float, float);
};

Atom::Atom()
{
   Symbol="AAA"; x[0]=0.0; x[1]=0.0; x[2]=0.0; R=0;
}

void Atom::read_Atom(string _Symbol, float _x, float _y, float _z)
{
   Symbol = _Symbol;
   x[0] = _x;
   x[1] = _y;
   x[2] = _z;
}

double Atomic_Distance(Atom atomo1, Atom atomo2)
{
   double suma=0;
   for(int ll=0;ll<3;ll++)
   {
      suma=suma+pow((atomo1.x[ll]-atomo2.x[ll]),2);
   }
   return sqrt(suma);
}

/********************************************************************/
/*********************** Molecule Definition ************************/
/********************************************************************/

class Molecule
{
   public:
      int Nat;
      Atom *atom;
      Molecule(string);
      Molecule();
      ~Molecule();
      void print_xyz(string);
      void move(float, float, float);
      void random_generator(string, int, string, int, float);
      void read_xyz(string);
      double x_min(); double x_max();
      double y_min(); double y_max();
      double z_min(); double z_max();
      bool fit_in(float, float, float, float, float, float);

};


/********************************************************************/
/********************** Emptiness Constructor ***********************/
/********************* Molecule  molecule_name; *********************/
/************************* Molecule Cysteine; ***********************/
/********************************************************************/

Molecule::Molecule()
{
  Nat=0;
  atom=NULL;
}

Molecule::~Molecule(){}

/********************************************************************/
/************************ Minimun_Separation ************************/
/******************** if(Minimum_Separation<2.0) ********************/
/********************************************************************/

double Minimun_Separation(Molecule Molecule1, Molecule Molecule2)
{
   double  min;
   int tot=Molecule1.Nat*Molecule2.Nat;
   int j;
   int k=0, l=0;
   double distances[tot];
   ofstream distancias("Distancias");

   for(i=0;i<Molecule1.Nat;i++)
   {
      for(j=0;j<Molecule2.Nat;j++)
      {
         distances[k]=Atomic_Distance(Molecule1.atom[i],Molecule2.atom[j])-Molecule1.atom[i].R-Molecule2.atom[j].R;
         distancias<<distances[k]<<endl;
         k++;
      }
   }
   distancias.close();
   system("cat Distancias | sort -n | head -1 >>Minimum_Separation");
   ifstream Minimum("Minimum_Separation");
   Minimum>>min;
   system("rm Minimum_Separation");
  return min;
}


/********************************************************************/
/******************************* z_min ******************************/
/************************* get_zmin(molecule) ***********************/
/********************************************************************/

double Molecule::z_max()
{
   double current, maximum=atom[0].x[2];
   for (i=1;i<Nat;i++)
   {
      current=atom[i].x[2];
      if( current > maximum )
      {
         maximum=current;
      }
   }
   return maximum;
}


double Molecule::y_max()
{
   double current, maximum=atom[0].x[1];
   for (i=1;i<Nat;i++)
   {
      current=atom[i].x[1];
      if( current > maximum )
      {
         maximum=current;
      }
   }
   return maximum;
}


double Molecule::x_max()
{
   double current, maximum=atom[0].x[0];
   for (i=1;i<Nat;i++)
   {
      current=atom[i].x[0];
      if( current > maximum )
      {
         maximum=current;
      }
   }
   return maximum;
}
double Molecule::z_min()
{
   double current, minimum=atom[0].x[2];
   for (i=1;i<Nat;i++)
   {
      current=atom[i].x[2];
      if( current < minimum )
      {
         minimum=current;
      }
   }
   return minimum;
}

double Molecule::y_min()
{
   double current, minimum=atom[0].x[1];
   for (i=1;i<Nat;i++)
   {
      current=atom[i].x[1];
      if( current < minimum )
      {
         minimum=current;
      }
   }
   return minimum;
}

double Molecule::x_min()
{
   double current, minimum=atom[0].x[0];
   for (i=1;i<Nat;i++)
   {
      current=atom[i].x[0];
      if( current < minimum )
      {
         minimum=current;
      }
   }
   return minimum;
}

bool Molecule::fit_in(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax)
{
   float xrange=xmax-xmin;
   float yrange=ymax-ymin;
   float zrange=zmax-zmin;
   bool fit;
   float mol_xrange=x_max()-x_min();
   float mol_yrange=y_max()-y_min();
   float mol_zrange=z_max()-z_min();

   if(mol_xrange < xrange && mol_yrange < yrange && mol_zrange < zrange)
   {
      fit=true;
   }
   else
   {
      fit=false;
   }
   return fit;
}

/********************************************************************/
/************************ Radii Constructor *************************/
/******************** map<string, double> Radios; *******************/
/********************* Radios=radii_dictionary(); *******************/
/********************************************************************/

map<string, double> radii_dictionary()
{
  typedef pair<string, double> radio_atomico;

     map<string, double> Radii;

     Radii.insert( radio_atomico("Au", 1.44) );
     Radii.insert( radio_atomico("Ag", 1.66) );
     Radii.insert( radio_atomico("Ir", 1.35) );
     Radii.insert( radio_atomico("Cu",1.28 ) );
     Radii.insert( radio_atomico("Rh",1.34 ) );
     Radii.insert( radio_atomico("Pt",1.39 ) );
     Radii.insert( radio_atomico("Ce",1.82 ) );
     Radii.insert( radio_atomico("H" ,0.46 ) );
     Radii.insert( radio_atomico("O" ,0.74 ) );
     Radii.insert( radio_atomico("N" ,0.74 ) );
     Radii.insert( radio_atomico("S" ,1.04 ) );
     Radii.insert( radio_atomico("C" ,0.77 ) );
     Radii.insert( radio_atomico("P" ,1.10 ) );

    return Radii;
}


/********************************************************************/
/************************* Radii Assignment *************************/
/********** atom[i].R=assign_radii(Radios,atom[i].Symbol); **********/
/********************************************************************/

double assign_radii(map<string, double> _Radios, string _Symbol)
{
     map<string, double>::iterator p = _Radios.find(_Symbol);
     double _radii;
     if(p != _Radios.end())
     _radii= p->second;
     return _radii;
}

/********************************************************************/
/****************** Constructor (from .xyz file) ********************/
/************** Molecule  molecule_name(string file); ***************/
/******************* Molecule Cysteine(argv[1]); ********************/
/********************************************************************/

Molecule::Molecule(string file)
{
   float x,y,z;
   string Symbol;
   string command="head -1 ";
          command+=file;
          command+=" >> Nat";
   system(command.c_str());

   ifstream Nat_file("Nat");
            Nat_file >> Nat;
            Nat_file.close();

   system("rm Nat");
          command.clear();
          command="head -";
          command+=patch::to_string(Nat+2);
          command+=" ";
          command+=file;
          command+=" | tail -";
          command+=patch::to_string(Nat);
          command+=" >> coordinatesAux ";
   system(command.c_str());

   ifstream coordinates_file("coordinatesAux");
   atom=new Atom[Nat+1];
   i=0;
   while(!coordinates_file.eof())
   {
      coordinates_file>>Symbol>>x>>y>>z;
      atom[i].read_Atom(Symbol,x,y,z);
      i++;
   }
   coordinates_file.close();
   system("rm coordinatesAux");
   map<string, double> Radios;
   Radios=radii_dictionary();
   for(i=0;i<Nat;i++)
   {
      atom[i].R=assign_radii(Radios,atom[i].Symbol);
   }

}


/***************************** read_xyz *****************************/
/********************* Molecule  molecule_name; *********************/
/***************** molecule_name.read_xyz(argv[1]); *****************/
/********************************************************************/

void Molecule::read_xyz(string file)
{
   float x,y,z;
   string Symbol;
   string command="head -1 ";
          command+=file;
          command+=" >> Nat";
   system(command.c_str());

   ifstream Nat_file("Nat");
            Nat_file >> Nat;
            Nat_file.close();

   system("rm Nat");
          command.clear();
          command="head -";
          command+=patch::to_string(Nat+2);
          command+=" ";
          command+=file;
          command+=" | tail -";
          command+=patch::to_string(Nat);
          command+=" >> coordinatesAux ";
   system(command.c_str());

   ifstream coordinates_file("coordinatesAux");
   atom=new Atom[Nat+1];
   i=0;
   while(!coordinates_file.eof())
   {
      coordinates_file>>Symbol>>x>>y>>z;
      atom[i].read_Atom(Symbol,x,y,z);
      i++;
   }
   coordinates_file.close();
   system("rm coordinatesAux");
   map<string, double> Radios;
   Radios=radii_dictionary();
   for(i=0;i<Nat;i++)
   {
      atom[i].R=assign_radii(Radios,atom[i].Symbol);
   }
}

/**************************** print_xyz *****************************/
/********************* Molecule  molecule_name; *********************/
/**************** molecule_name.print_xyz(argv[1]); *****************/
/********************************************************************/

void Molecule::print_xyz(string outputfile)
{
   FILE *p = fopen(outputfile.c_str(),"w");
   fprintf(p,"%i \n \n", Nat);
   for(i=0;i<Nat;i++)
   {
      fprintf(p,"%s  %f  %f  %f \n", atom[i].Symbol.c_str(),atom[i].x[0],atom[i].x[1],atom[i].x[2]);
   }
fclose(p);
/*
   ofstream output_file(outputfile);
   output_file<<Nat<<endl;
   output_file<<" "<<endl;
   for(i=0;i<Nat;i++)
   {
      output_file<<atom[i].Symbol<<" "<<atom[i].x[0]<<" "<<atom[i].x[1]<<" "<<atom[i].x[2]<<endl;
   }
   output_file.close();
*/
}

/************************* random_generator *************************/
/********************** Cluster  Cluster_name; **********************/
/******** Cluster_name.random_generator("Au",5,"Ir",3,range) ********/
/************* Cluster_name.random_generator("Ir",3) ****************/
/********** Cluster_name.random_generator("Ir",3,range) *************/


void Molecule::random_generator(string Symbol_1, int N_Symbol_1, string Symbol_2="AAA", int N_Symbol_2=0, float epsilon=2.5)
{

  map<string, double> Radios;
   Radios=radii_dictionary();
   Nat=N_Symbol_1+N_Symbol_2;
   int random;
   int randomS;
   double criterio;
   int accepted=0;
   int rejected=0;
   float Mx,Nx,My,Ny,Mz,Nz;
   float x,y,z;
   double Distance;
   srand(time(NULL));
   atom=new Atom[Nat+1];
   int cont_S1=1;
   int cont_S2=0;

///////////////////////////////Coloca el primer
///////////////////////////////átomo en el origen
   atom[0].Symbol=Symbol_1;  //o podría ser que
   atom[0].x[0]=0;           //se genere en un
   atom[0].x[1]=0;           //punto aleatorio
   atom[0].x[2]=0;           //o que admita como
///////////////////////////////argumento este punto
///////////////////////////////por default el origen


   atom[0].R=assign_radii(Radios, atom[0].Symbol);
   for(i=1;i<Nat;i++)
   {
      accepted=0;
      while(accepted==0)
      {
         random=(i-1)+(double)rand()/((double)RAND_MAX/(0-(i-1)+1)+1);

         if(N_Symbol_2!=0)
         {
            randomS=rand()%2;

            if(cont_S1<N_Symbol_1 && randomS==0)
            {
               atom[i].Symbol=Symbol_1;
            }
            else
            {
               if(cont_S2<N_Symbol_2 && randomS==1)
               {
                  atom[i].Symbol=Symbol_2;
               }
               else
               {
                  if(cont_S1>=N_Symbol_1)
                  {
                     atom[i].Symbol=Symbol_2;
                  }
                  else
                  {
                     if(cont_S2>=N_Symbol_2)
                     {
                        atom[i].Symbol=Symbol_1;
                     }
                  }

               }
            }
         }
         else
         {
            atom[i].Symbol=Symbol_1;
         }
         atom[i].R=assign_radii(Radios, atom[i].Symbol);

         Mx=atom[random].x[0]+epsilon;  //Maximum range for random X
         Nx=atom[random].x[0]-epsilon;  //Minimun range for random x
         My=atom[random].x[1]+epsilon;  //Maximum range for random y
         Ny=atom[random].x[1]-epsilon;  //Minimum range for random y
         Mz=atom[random].x[2]+epsilon;  //Maximum range for random z
         Nz=atom[random].x[2]-epsilon;  //Minimum range for randon z

         atom[i].x[0]=Nx+ ((double)rand())/((double)RAND_MAX )* (Mx-Nx);
         atom[i].x[1]=Ny+ ((double)rand())/((double)RAND_MAX )* (My-Ny);
         atom[i].x[2]=Nz+ ((double)rand())/((double)RAND_MAX )* (Mz-Nz);

         rejected=0;

         for(j=0;j<i;j++)
         {
            Distance=sqrt(pow(atom[i].x[0]-atom[j].x[0],2)+pow(atom[i].x[1]-atom[j].x[1],2)+pow(atom[i].x[2]-atom[j].x[2],2));
            criterio=atom[i].R+atom[j].R;
            if(Distance<criterio)
            {
               rejected++;
            }
         }

         if(rejected>0)
         {
            accepted=0;
         }
         else
         {
            accepted=1;
            if(strcmp(atom[i].Symbol.c_str(),Symbol_1.c_str()) == 0 )
            {
               cont_S1++;
            }
            else
            {
               cont_S2++;
            }

         }
         if(i<Nat)
         {
            continue;
         }
         else
         {
            break;
         }

      }

   }

}


/*****************************  move ********************************/
/************** molecule_name.move(DeltaX,DeltaY,DeltaZ) ************/
/********************** Cysteine.move(1,5,3) ************************/
/********************************************************************/

void Molecule::move(float Dx, float Dy, float Dz)
{
   for(i=0;i<Nat;i++)
   {
       atom[i].x[0]=atom[i].x[0]+Dx;
       atom[i].x[1]=atom[i].x[1]+Dy;
       atom[i].x[2]=atom[i].x[2]+Dz;
   }
}

float random_number(float min,float max)
{
   float random = min + ((float)rand())/((float)RAND_MAX ) * (max-min);
   return random;
}

int main(int argc, char *argv[])
{
double _x, _y, _z;
Molecule clus;
if(argc==11)
{
  clus.random_generator(argv[1],atoi(argv[2]));
  bool acep=clus.fit_in(atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[9]),atoi(argv[10]));
  while (acep==false)
  {
     clus.~Molecule();
     clus.random_generator(argv[1],atoi(argv[2]));
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

clus.random_generator(argv[1],atoi(argv[2]),argv[3],atoi(argv[4]));
bool acep=clus.fit_in(atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),atoi(argv[11]),atoi(argv[12]));
while (acep==false)
{
   clus.~Molecule();
   clus.random_generator(argv[1],atoi(argv[2]),argv[3],atoi(argv[4]));
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
