//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/ atomicpp.h library _/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_ For simulation of atoms and molecules _/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
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
#include   <cmath>
#include   <ctime>


#include <string>
#include <utility>
#include <map>
#include <iomanip>

/********* * * * * * *  *  *   *   *   *  *  *  * * * ***************/

int Nat;
int i,j,k,ll;
using namespace std;


/********************************************************************/
/************************* Atom Definition **************************/
/********************************************************************/

class Atom
{
   public:
     string Symbol;
     double x[3];
     double v[3];
     double a[3];
     double R;
     Atom();
     void read_Atom(string, float, float, float);
};

Atom::Atom()
{
   Symbol="AAA"; R=0;
   x[0]=0.0; x[1]=0.0; x[2]=0.0;
   v[0]=0.0; v[1]=0.0; v[2]=0.0;
   a[0]=0.0; a[1]=0.0; a[2]=0.0;
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
   for(ll=0;ll<3;ll++)
   {
      suma=suma+pow(atomo1.x[ll]-atomo2.x[ll],2);
   }
   return sqrt(suma);
}


/********************************************************************/
/******************* Atomic_Structure Definition ********************/
/********************************************************************/

class Atomic_Structure{
     public:
        int Nat;
        Atom *atom;
        Atomic_Structure(string);
        Atomic_Structure();
        ~Atomic_Structure();
        void read_xyz(string file);
        void read_fhi(string file);
        void read_VASP(string file);
        void print_xyz(string);
        void print_VASP(string, string, float, float (*)[3]);
        double x_min(); double x_max();
        double y_min(); double y_max();
        double z_min(); double z_max();
        bool fit_in(float, float, float, float, float, float);
};

/********************************************************************/
/************************ Cluster Definition ************************/
/********************************************************************/


class Cluster : public Atomic_Structure{
     public:
        string type;
        double quirality;
        Cluster(string);
        Cluster();
        ~Cluster();
        void move(float, float, float);
        void rotate_Rad(float, float);
        void rotate_Deg(float, float);
        void kick(float);
        void kick_lennard(float);
        void swap(int);
        void srand_generator(string, int, string, int, float);
        void rand_generator(string, int, string, int);
        void centroid();
};


/********************************************************************/
/************************ Crystal Definition ************************/
/********************************************************************/


class Crystal : public Atomic_Structure{
     public:
       //Parámetros de celda
       double x[3][3];
       double factor;
        Crystal(string);
        Crystal();
        ~Crystal();
};

/********************************************************************/
/*********************** Molecule Definition ************************/
/********************************************************************/

class Molecule : public Atomic_Structure{
    public:
       Molecule(string);
       Molecule();
       ~Molecule();
      void move(float, float, float);
      void rotate_Rad(float, float);
      void rotate_Deg(float, float);
      void centroid();

};


/**************************** random_number *************************/
/********************* float r=random_number(0,10); *****************/
/********************************************************************/

float random_number(float min,float max)
{
   float random = min + ((float)rand())/((float)RAND_MAX ) * (max-min);
   return random;
}


/********************************************************************/
/********************** Emptiness Constructor ***********************/
/***************** Atomic_Structure  molecule_name; *****************/
/************************* Molecule Cysteine; ***********************/
/********************************************************************/

Atomic_Structure::Atomic_Structure()
{
  Nat=0;
  atom=NULL;
}

/********************************************************************/
/********************** Emptiness Constructor ***********************/
/********************* Molecule  molecule_name; *********************/
/************************* Molecule Cysteine; ***********************/
/********************************************************************/

Molecule::Molecule() : Atomic_Structure()
{
  Nat=0;
  atom=NULL;
}

/********************************************************************/
/********************** Emptiness Constructor ***********************/
/********************** Cluster  cluster_name; **********************/
/************************* Molecule Cysteine; ***********************/
/********************************************************************/

Cluster::Cluster() : Atomic_Structure()
{
  Nat=0;
  atom=NULL;
}

/********************************************************************/
/********************** Emptiness Constructor ***********************/
/********************** Crystal  crystal_name; **********************/
/************************* Molecule Cysteine; ***********************/
/********************************************************************/

Crystal::Crystal() : Atomic_Structure()
{
  Nat=0;
  atom=NULL;
}

Atomic_Structure::~Atomic_Structure(){}
Molecule::~Molecule(){}
Cluster::~Cluster(){}
Crystal::~Crystal(){}

/********************************************************************/
/************************ minimum_separation ************************/
/**************** Molecule  molecule_1, molecule_2; *****************/
/******* separation=minimum_separation(molecule_1,molecule_2); ******/
/********************************************************************/

double minimun_separation(Atomic_Structure Molecule1, Atomic_Structure Molecule2)
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

Atomic_Structure::Atomic_Structure(string file)
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
          command+=to_string(Nat+2);
          command+=" ";
          command+=file;
          command+=" | tail -";
          command+=to_string(Nat);
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

/********************************************************************/
/******************************* z_max ******************************/
/************************** molecule.z_max() ************************/
/********************************************************************/

double Atomic_Structure::z_max()
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


double Atomic_Structure::y_max()
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


double Atomic_Structure::x_max()
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
double Atomic_Structure::z_min()
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

double Atomic_Structure::y_min()
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

double Atomic_Structure::x_min()
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


/******************************* fit_in *****************************/
/*********** molecule.fit_in(xmin,xmax,ymin,ymax,zmin,zmax **********/
/********************************************************************/

bool Atomic_Structure::fit_in(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax)
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


/***************************** read_xyz *****************************/
/********************* Molecule  molecule_name; *********************/
/***************** molecule_name.read_xyz(argv[1]); *****************/
/********************************************************************/

void Atomic_Structure::read_xyz(string file)
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
          command+=to_string(Nat+2);
          command+=" ";
          command+=file;
          command+=" | tail -";
          command+=to_string(Nat);
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

/***************************** read_fhi *****************************/
/********************* Molecule  molecule_name; *********************/
/**************** molecule_name.read_fhi("file.in"); ****************/
/********************************************************************/

void Atomic_Structure::read_fhi(string file)
{
   string com="cp ";
          com+=file;
          com+=" tmp.in";
   system(com.c_str());
   system("grep \"atom\" tmp.in | wc -l > tmp.xyz");
   system("echo \" \" >> tmp.xyz");
   system("grep \"atom\" geometry.in | awk '{print $5\" \"$2\" \"$3\" \"$4}' >> tmp.xyz");
   system("rm tmp.in");
   float x,y,z;
   string Symbol;

   string command="head -1 ";
          command+="tmp.xyz";
          command+=" >> Nat";
   system(command.c_str());

   ifstream Nat_file("Nat");
            Nat_file >> Nat;
            Nat_file.close();

   system("rm Nat");
          command.clear();
          command="head -";
          command+=to_string(Nat+2);
          command+=" ";
          command+="tmp.xyz";
          command+=" | tail -";
          command+=to_string(Nat);
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


/*************************** rand_generator *************************/
/********************** Cluster  Cluster_name; **********************/
/********** Cluster_name.rand_generator("Au",5,"Ir",3,range) ********/
/*************** Cluster_name.rand_generator("Ir",3) ****************/
/************ Cluster_name.rand_generator("Ir",3,range) *************/


void Cluster::rand_generator(string Symbol_1, int N_Symbol_1, string Symbol_2="AAA", int N_Symbol_2=0 )
{
   map<string, double> Radios;
   Radios=radii_dictionary();
   Nat=N_Symbol_1+N_Symbol_2;
   srand(time(NULL));
   atom=new Atom[Nat+1];
   int cont_S1=1;
   int cont_S2=0;
   int randomS;

   if(strcmp(Symbol_2.c_str(),"AAA") != 0 )
   {
      type="bimetallic";
   }
   else
   {
      type="monometallic";
   }
   double pseudo_rand[81][3];
   pseudo_rand[1][0]=-0.0195029 ; pseudo_rand[1][1]=0.0847575 ; pseudo_rand[1][2]=-0.0264053;
   pseudo_rand[2][0]=-0.029306 ; pseudo_rand[2][1]=-2.51719 ; pseudo_rand[2][2]=-0.0385384;
   pseudo_rand[3][0]=1.32679 ; pseudo_rand[3][1]=-1.18825 ; pseudo_rand[3][2]=1.79442;
   pseudo_rand[4][0]=1.24944 ; pseudo_rand[4][1]=-1.2225 ; pseudo_rand[4][2]=-1.89543;
   pseudo_rand[5][0]=-2.17871 ; pseudo_rand[5][1]=-1.22141 ; pseudo_rand[5][2]=0.640053;
   pseudo_rand[7][0]=-0.838767 ; pseudo_rand[7][1]=0.0835188 ; pseudo_rand[7][2]=2.45218;
   pseudo_rand[8][0]=2.59034 ; pseudo_rand[8][1]=0.0824258 ; pseudo_rand[8][2]=-0.0840119;
   pseudo_rand[6][0]=-1.35706 ; pseudo_rand[6][1]=-1.23125 ; pseudo_rand[6][2]=-1.83431;
   pseudo_rand[9][0]=1.32544 ; pseudo_rand[9][1]=1.38315 ; pseudo_rand[9][2]=1.79148;
   pseudo_rand[10][0]=1.24746 ; pseudo_rand[10][1]=1.38283 ; pseudo_rand[10][2]=-1.89791;
   pseudo_rand[11][0]=-2.17943 ; pseudo_rand[11][1]=1.38392 ; pseudo_rand[11][2]=0.636637;
   pseudo_rand[12][0]=-1.35706 ; pseudo_rand[12][1]=1.39101 ; pseudo_rand[12][2]=-1.83544;
   pseudo_rand[13][0]=-0.0259896 ; pseudo_rand[13][1]=2.68706 ; pseudo_rand[13][2]=-0.0362973;
   pseudo_rand[14][0]=-0.829979 ; pseudo_rand[14][1]=-2.48665 ; pseudo_rand[14][2]=2.42262;
   pseudo_rand[15][0]=2.55844 ; pseudo_rand[15][1]=-2.48773 ; pseudo_rand[15][2]=-0.0834814;
   pseudo_rand[16][0]=-0.0928186 ; pseudo_rand[16][1]=0.077041 ; pseudo_rand[16][2]=-3.65928;
   pseudo_rand[17][0]=-3.47169 ; pseudo_rand[17][1]=0.078118 ; pseudo_rand[17][2]=-1.16024;
   pseudo_rand[18][0]=-0.826892 ; pseudo_rand[18][1]=2.65276 ; pseudo_rand[18][2]=2.41527;
   pseudo_rand[19][0]=2.55262 ; pseudo_rand[19][1]=2.65168 ; pseudo_rand[19][2]=-0.0842503;
   pseudo_rand[20][0]=1.31501 ; pseudo_rand[20][1]=-3.74875 ; pseudo_rand[20][2]=1.7796;
   pseudo_rand[22][0]=-2.16839 ; pseudo_rand[22][1]=-3.82939 ; pseudo_rand[22][2]=0.623685;
   pseudo_rand[23][0]=-1.356 ; pseudo_rand[23][1]=-3.82677 ; pseudo_rand[23][2]=-1.83176;
   pseudo_rand[21][0]=1.22969 ; pseudo_rand[21][1]=-3.83047 ; pseudo_rand[21][2]=-1.88957;
   pseudo_rand[24][0]=0.524113 ; pseudo_rand[24][1]=-1.2125 ; pseudo_rand[24][2]=4.26029;
   pseudo_rand[25][0]=3.9196 ; pseudo_rand[25][1]=-1.21359 ; pseudo_rand[25][2]=1.74896;
   pseudo_rand[26][0]=-0.099016 ; pseudo_rand[26][1]=-2.54239 ; pseudo_rand[26][2]=-3.68792;
   pseudo_rand[27][0]=-3.49834 ; pseudo_rand[27][1]=-2.54131 ; pseudo_rand[27][2]=-1.17375;
   pseudo_rand[28][0]=2.66405 ; pseudo_rand[28][1]=0.0885545 ; pseudo_rand[28][2]=3.60194;
   pseudo_rand[29][0]=-2.99977 ; pseudo_rand[29][1]=-1.21263 ; pseudo_rand[29][2]=3.1135;
   pseudo_rand[30][0]=3.8549 ; pseudo_rand[30][1]=-1.21482 ; pseudo_rand[30][2]=-1.95627;
   pseudo_rand[31][0]=0.521314 ; pseudo_rand[31][1]=1.36949 ; pseudo_rand[31][2]=4.25938;
   pseudo_rand[32][0]=3.92062 ; pseudo_rand[32][1]=1.3684 ; pseudo_rand[32][2]=1.74523;
   pseudo_rand[33][0]=2.52371 ; pseudo_rand[33][1]=0.0752852 ; pseudo_rand[33][2]=-3.75257;
   pseudo_rand[34][0]=-4.32703 ; pseudo_rand[34][1]=0.0774689 ; pseudo_rand[34][2]=1.3143;
   pseudo_rand[36][0]=3.85082 ; pseudo_rand[36][1]=1.36547 ; pseudo_rand[36][2]=-1.95956;
   pseudo_rand[35][0]=-3.00066 ; pseudo_rand[35][1]=1.36766 ; pseudo_rand[35][2]=3.10785;
   pseudo_rand[37][0]=-2.69634 ; pseudo_rand[37][1]=0.0758621 ; pseudo_rand[37][2]=-3.64567;
   pseudo_rand[39][0]=-0.101921 ; pseudo_rand[39][1]=2.6949 ; pseudo_rand[39][2]=-3.69677;
   pseudo_rand[38][0]=1.31946 ; pseudo_rand[38][1]=3.98582 ; pseudo_rand[38][2]=1.78228;
   pseudo_rand[40][0]=-3.5038 ; pseudo_rand[40][1]=2.69599 ; pseudo_rand[40][2]=-1.18072;
   pseudo_rand[41][0]=-2.17698 ; pseudo_rand[41][1]=3.98769 ; pseudo_rand[41][2]=0.613141;
   pseudo_rand[42][0]=1.22535 ; pseudo_rand[42][1]=3.98661 ; pseudo_rand[42][2]=-1.90325;
   pseudo_rand[43][0]=-1.35967 ; pseudo_rand[43][1]=3.98151 ; pseudo_rand[43][2]=-1.84008;
   pseudo_rand[44][0]=-0.00349797 ; pseudo_rand[44][1]=-5.09938 ; pseudo_rand[44][2]=-0.00253142;
   pseudo_rand[45][0]=2.64629 ; pseudo_rand[45][1]=-2.50334 ; pseudo_rand[45][2]=3.57904;
   pseudo_rand[46][0]=2.48502 ; pseudo_rand[46][1]=-2.51022 ; pseudo_rand[46][2]=-3.72354;
   pseudo_rand[47][0]=-4.28901 ; pseudo_rand[47][1]=-2.50806 ; pseudo_rand[47][2]=1.28659;
   pseudo_rand[48][0]=-1.62841 ; pseudo_rand[48][1]=0.0750585 ; pseudo_rand[48][2]=4.88148;
   pseudo_rand[49][0]=5.14437 ; pseudo_rand[49][1]=0.0728996 ; pseudo_rand[49][2]=-0.127727;
   pseudo_rand[50][0]=-2.68054 ; pseudo_rand[50][1]=-2.49415 ; pseudo_rand[50][2]=-3.62319;
   pseudo_rand[51][0]=2.63614 ; pseudo_rand[51][1]=2.68213 ; pseudo_rand[51][2]=3.56308;
   pseudo_rand[52][0]=-4.29393 ; pseudo_rand[52][1]=2.66266 ; pseudo_rand[52][2]=1.28107;
   pseudo_rand[53][0]=2.48332 ; pseudo_rand[53][1]=2.6605 ; pseudo_rand[53][2]=-3.73144;
   pseudo_rand[54][0]=-2.68333 ; pseudo_rand[54][1]=2.6446 ; pseudo_rand[54][2]=-3.62919;
   pseudo_rand[55][0]=-0.0111239 ; pseudo_rand[55][1]=5.24589 ; pseudo_rand[55][2]=-0.0173006;
   pseudo_rand[57][0]=3.86127 ; pseudo_rand[57][1]=-3.78637 ; pseudo_rand[57][2]=1.71154;
   pseudo_rand[58][0]=-2.93427 ; pseudo_rand[58][1]=-3.76959 ; pseudo_rand[58][2]=3.06773;
   pseudo_rand[56][0]=0.504357 ; pseudo_rand[56][1]=-3.7853 ; pseudo_rand[56][2]=4.19434;
   pseudo_rand[59][0]=3.7909 ; pseudo_rand[59][1]=-3.77173 ; pseudo_rand[59][2]=-1.90626;
   pseudo_rand[60][0]=-0.82015 ; pseudo_rand[60][1]=-5.06436 ; pseudo_rand[60][2]=2.43582;
   pseudo_rand[62][0]=-1.61014 ; pseudo_rand[62][1]=-2.4979 ; pseudo_rand[62][2]=4.84251;
   pseudo_rand[61][0]=2.56712 ; pseudo_rand[61][1]=-5.06544 ; pseudo_rand[61][2]=-0.0694326;
   pseudo_rand[64][0]=1.14477 ; pseudo_rand[64][1]=-1.20367 ; pseudo_rand[64][2]=-5.4834;
   pseudo_rand[63][0]=5.1007 ; pseudo_rand[63][1]=-2.50004 ; pseudo_rand[63][2]=-0.120882;
   pseudo_rand[65][0]=-5.57876 ; pseudo_rand[65][1]=-1.20153 ; pseudo_rand[65][2]=-0.510616;
   pseudo_rand[66][0]=-1.40233 ; pseudo_rand[66][1]=-1.20338 ; pseudo_rand[66][2]=-5.44063;
   pseudo_rand[68][0]=-1.62222 ; pseudo_rand[68][1]=2.64533 ; pseudo_rand[68][2]=4.83796;
   pseudo_rand[69][0]=5.10201 ; pseudo_rand[69][1]=2.64319 ; pseudo_rand[69][2]=-0.135334;
   pseudo_rand[70][0]=1.14364 ; pseudo_rand[70][1]=1.35008 ; pseudo_rand[70][2]=-5.48765;
   pseudo_rand[71][0]=-5.58144 ; pseudo_rand[71][1]=1.35222 ; pseudo_rand[71][2]=-0.513724;
   pseudo_rand[67][0]=-4.79206 ; pseudo_rand[67][1]=-1.2023 ; pseudo_rand[67][2]=-2.93356;
   pseudo_rand[72][0]=3.87373 ; pseudo_rand[72][1]=3.92977 ; pseudo_rand[72][2]=1.68879;
   pseudo_rand[73][0]=0.482131 ; pseudo_rand[73][1]=3.93085 ; pseudo_rand[73][2]=4.19725;
   pseudo_rand[74][0]=-2.93769 ; pseudo_rand[74][1]=3.92401 ; pseudo_rand[74][2]=3.06034;
   pseudo_rand[75][0]=3.78801 ; pseudo_rand[75][1]=3.92187 ; pseudo_rand[75][2]=-1.91404;
   pseudo_rand[76][0]=-1.40406 ; pseudo_rand[76][1]=1.35043 ; pseudo_rand[76][2]=-5.44481;
   pseudo_rand[77][0]=-4.79449 ; pseudo_rand[77][1]=1.35151 ; pseudo_rand[77][2]=-2.93722;
   pseudo_rand[78][0]=-0.832405 ; pseudo_rand[78][1]=5.20814 ; pseudo_rand[78][2]=2.41864;
   pseudo_rand[79][0]=2.55852 ; pseudo_rand[79][1]=5.20706 ; pseudo_rand[79][2]=-0.0893115;
   pseudo_rand[80][0]=1.36354 ; pseudo_rand[80][1]=-6.22286 ; pseudo_rand[80][2]=1.84628;

   atom[0].Symbol=Symbol_1;  //o podría ser que
   atom[0].x[0]=pseudo_rand[1][0];           //se genere en un
   atom[0].x[1]=pseudo_rand[1][1];           //punto aleatorio
   atom[0].x[2]=pseudo_rand[1][2];
   atom[0].R=assign_radii(Radios, atom[0].Symbol);

   for(i=1;i<Nat;i++)
   {
      if(N_Symbol_2!=0)
      {
         randomS=rand()%2;
         if(cont_S1<N_Symbol_1 && randomS==0)
         {
            atom[i].Symbol=Symbol_1;
            cont_S1++;
         }
         else
         {
            if(cont_S2<N_Symbol_2 && randomS==1)
            {
               atom[i].Symbol=Symbol_2;
               cont_S2++;
            }
            else
            {
               if(cont_S1>=N_Symbol_1)
               {
                  atom[i].Symbol=Symbol_2;
                  cont_S2++;
               }
               else
               {
                  if(cont_S2>=N_Symbol_2)
                  {
                     atom[i].Symbol=Symbol_1;
                     cont_S1++;
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

      for(j=0;j<3;j++)
      {
         atom[i].x[j]=pseudo_rand[i+1][j];
      }
   }
}



/************************* srand_generator **************************/
/********************** Cluster  Cluster_name; **********************/
/********* Cluster_name.srand_generator("Au",5,"Ir",3,range) ********/
/************** Cluster_name.srand_generator("Ir",3) ****************/
/*********** Cluster_name.srand_generator("Ir",3,range) *************/


void Cluster::srand_generator(string Symbol_1, int N_Symbol_1, string Symbol_2="AAA", int N_Symbol_2=0, float epsilon=2.5)
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

   if(strcmp(Symbol_2.c_str(),"AAA") != 0 )
   {
      type="bimetallic";
   }
   else
   {
      type="monometallic";
   }
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



/************************** print__xyz ******************************/
/********************* Molecule  molecule_name; *********************/
/**************** molecule_name.print_xyz(argv[1]); *****************/
/********************************************************************/

void Atomic_Structure::print_xyz(string outputfile)
{
   ofstream output_file(outputfile);
   output_file<<Nat<<endl;
   output_file<<" "<<endl;
   for(i=0;i<Nat;i++)
   {
      output_file<<atom[i].Symbol<<" "<<atom[i].x[0]<<" "<<atom[i].x[1]<<" "<<atom[i].x[2]<<endl;
   }
   output_file.close();
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

/*****************************  move ********************************/
/************** molecule_name.move(DeltaX,DeltaY,DeltaZ) ************/
/********************** Cysteine.move(1,5,3) ************************/
/********************************************************************/

void Cluster::move(float Dx, float Dy, float Dz)
{
   for(i=0;i<Nat;i++)
   {
       atom[i].x[0]=atom[i].x[0]+Dx;
       atom[i].x[1]=atom[i].x[1]+Dy;
       atom[i].x[2]=atom[i].x[2]+Dz;
   }
}

/*************************** rotate_Rad *****************************/
/********** molecule_name.rotate_Rad(theta(rad), phi(rad)) **********/
/****************** Cysteine.rotate_Rad(pi,pi/2) ********************/
/********************************************************************/

void Molecule::rotate_Rad(float theta, float phi)
{
   float s;
   int l,m;
   float Rot[3][3];

   Rot[0][0]=cos(phi);                   //*************************//
   Rot[0][1]=(-cos(theta)*sin(phi));     //                         //
   Rot[0][2]=sin(phi)*sin(theta);        //                         //
   Rot[1][0]=sin(phi);                   //                         //
   Rot[1][1]=cos(theta)*cos(phi);        //     Rotation Matriz     //
   Rot[1][2]=(-sin(theta)*cos(phi));     //                         //
   Rot[2][0]=0;                          //                         //
   Rot[2][1]=sin(theta);                 //                         //
   Rot[2][2]=cos(theta);                 //*************************//

   for(i=0;i<Nat;i++)
   {
      for(l=0;l<3;l++)
      {
         s=0;
         for(m=0;m<3;m++)
         {
            s=s+(Rot[l][m]*(atom[i].x[m]));
         }

         atom[i].x[l]=s;
      }

   }

}

/*************************** rotate_Rad *****************************/
/********** molecule_name.rotate_Rad(theta(rad), phi(rad)) **********/
/****************** Cysteine.rotate_Rad(pi,pi/2) ********************/
/********************************************************************/

void Cluster::rotate_Rad(float theta, float phi)
{
   float s;
   int l,m;
   float Rot[3][3];

   Rot[0][0]=cos(phi);                   //*************************//
   Rot[0][1]=(-cos(theta)*sin(phi));     //                         //
   Rot[0][2]=sin(phi)*sin(theta);        //                         //
   Rot[1][0]=sin(phi);                   //                         //
   Rot[1][1]=cos(theta)*cos(phi);        //     Rotation Matriz     //
   Rot[1][2]=(-sin(theta)*cos(phi));     //                         //
   Rot[2][0]=0;                          //                         //
   Rot[2][1]=sin(theta);                 //                         //
   Rot[2][2]=cos(theta);                 //*************************//

   for(i=0;i<Nat;i++)
   {
      for(l=0;l<3;l++)
      {
         s=0;
         for(m=0;m<3;m++)
         {
            s=s+(Rot[l][m]*(atom[i].x[m]));
         }

         atom[i].x[l]=s;
      }

   }

}

/*************************** rotate_Deg *****************************/
/********** molecule_name.rotate_Deg(theta(deg), phi(deg)) **********/
/******************* Cysteine.rotate_Deg(45,90) *********************/
/********************************************************************/

void Molecule::rotate_Deg(float thetadeg, float phideg)
{
   float s;
   int l,m;
   float Rot[3][3];
   float theta,phi;

   theta=(thetadeg*(3.1415926535/180));
   phi=(phideg*(3.1415926535/180));

   Rot[0][0]=cos(phi);                   //*************************//
   Rot[0][1]=(-cos(theta)*sin(phi));     //                         //
   Rot[0][2]=sin(phi)*sin(theta);        //                         //
   Rot[1][0]=sin(phi);                   //                         //
   Rot[1][1]=cos(theta)*cos(phi);        //     Rotation Matriz     //
   Rot[1][2]=(-sin(theta)*cos(phi));     //                         //
   Rot[2][0]=0;                          //                         //
   Rot[2][1]=sin(theta);                 //                         //
   Rot[2][2]=cos(theta);                 //*************************//

   for(i=0;i<Nat;i++)
   {
      for(l=0;l<3;l++)
      {
         s=0;
         for(m=0;m<3;m++)
         {
            s=s+(Rot[l][m]*(atom[i].x[m]));
         }
         atom[i].x[l]=s;
      }
   }

}


/*************************** rotate_Deg *****************************/
/********** molecule_name.rotate_Deg(theta(deg), phi(deg)) **********/
/******************* Cysteine.rotate_Deg(45,90) *********************/
/********************************************************************/

void Cluster::rotate_Deg(float thetadeg, float phideg)
{
   float s;
   int l,m;
   float Rot[3][3];
   float theta,phi;

   theta=(thetadeg*(3.1415926535/180));
   phi=(phideg*(3.1415926535/180));

   Rot[0][0]=cos(phi);                   //*************************//
   Rot[0][1]=(-cos(theta)*sin(phi));     //                         //
   Rot[0][2]=sin(phi)*sin(theta);        //                         //
   Rot[1][0]=sin(phi);                   //                         //
   Rot[1][1]=cos(theta)*cos(phi);        //     Rotation Matriz     //
   Rot[1][2]=(-sin(theta)*cos(phi));     //                         //
   Rot[2][0]=0;                          //                         //
   Rot[2][1]=sin(theta);                 //                         //
   Rot[2][2]=cos(theta);                 //*************************//

   for(i=0;i<Nat;i++)
   {
      for(l=0;l<3;l++)
      {
         s=0;
         for(m=0;m<3;m++)
         {
            s=s+(Rot[l][m]*(atom[i].x[m]));
         }
         atom[i].x[l]=s;
      }
   }

}

/******************************  Kick *******************************/
/***************** molecule_name.kick(step_width) *******************/
/*********************** Cysteine.kick(0.8) *************************/
/********************************************************************/

void Cluster::kick(float step_width)
{

   float DelX, DelY, DelZ;
   srand(time(NULL));

   for(i=0;i<Nat;i++)
   {
      DelX=((double)rand()/((double)RAND_MAX))*step_width;
      DelY=((double)rand()/((double)RAND_MAX))*step_width;
      DelZ=((double)rand()/((double)RAND_MAX))*step_width;

      atom[i].x[0]=DelX+atom[i].x[0];
      atom[i].x[1]=DelY+atom[i].x[1];
      atom[i].x[2]=DelZ+atom[i].x[2];
   }

}

/************************  Kick_Lennard *****************************/
/************ cluster_name.kick_lennard(step_width) *****************/
/******************* Cysteine.kick_lennard(0.8) *********************/
/********************************************************************/

void Cluster::kick_lennard(float kick_size=1.5)
{

   int time_criterio=0;
   float time=0;
   float time_step=0.0001;
   float G_=-1.0;
   //float kick_size=1.50;  //1.5 funciona bien para 40
   float ax,ay,az;
   for(i=0;i<Nat;i++)
   {
      atom[i].x[0]=atom[i].x[0]*kick_size+random_number(-1,1)*kick_size/2;
      atom[i].x[1]=atom[i].x[1]*kick_size+random_number(-1,1)*kick_size/2;
      atom[i].x[2]=atom[i].x[2]*kick_size+random_number(-1,1)*kick_size/2;
   }
   float Dist;
   float sep;
   double swap;
   while(time_criterio==0)
   {
      for(i=0;i<Nat;i++) //CALCULA LAS FUERZAS
      {
         ax=0;ay=0;az=0;
         for(j=0;j<Nat;j++)
         {
            if(i!=j)
            {
               //Calcula
               Dist=Atomic_Distance(atom[i],atom[j]);
               ax=ax+G_*(4*((12*pow((atom[i].R+atom[j].R)*0.95,12)/pow(Dist,13))-(6*pow((atom[i].R+atom[j].R)*0.95,6)/pow(Dist,7)))*(atom[i].x[0]-atom[j].x[0]));
               ay=ay+G_*(4*((12*pow((atom[i].R+atom[j].R)*0.95,12)/pow(Dist,13))-(6*pow((atom[i].R+atom[j].R)*0.95,6)/pow(Dist,7)))*(atom[i].x[1]-atom[j].x[1]));
               az=az+G_*(4*((12*pow((atom[i].R+atom[j].R)*0.95,12)/pow(Dist,13))-(6*pow((atom[i].R+atom[j].R)*0.95,6)/pow(Dist,7)))*(atom[i].x[2]-atom[j].x[2]));
            }
         }
         atom[i].a[0]=ax*G_;
         atom[i].a[1]=ay*G_;
         atom[i].a[2]=az*G_;
      }
      for(i=0;i<Nat;i++)
      {
         //// Velocidades
         atom[i].v[0]=atom[i].v[0]+(time_step*atom[i].a[0]);//time;
         atom[i].v[1]=atom[i].v[1]+(time_step*atom[i].a[1]);//time;
         atom[i].v[2]=atom[i].v[2]+(time_step*atom[i].a[2]);//time;
         //// Posiciones
         atom[i].x[0]=atom[i].x[0]+(time_step*atom[i].v[0]);
         atom[i].x[1]=atom[i].x[1]+(time_step*atom[i].v[1]);
         atom[i].x[2]=atom[i].x[2]+(time_step*atom[i].v[2]);
      }
      time++;
      if(time>12500)
      {
         time_criterio=1;
      }
   }
}

/******************************* swap *******************************/
/************************ swap(TimesToSwap) *************************/
/************************** Cluster.swap(3) *************************/

void Cluster::swap(int Nswap)
{
   int a,b,i;

   string aux, aux2;
   int M=0;
   int N=Nat-1;
   srand(time(NULL));


   i=0;

   if(strcmp(type.c_str(),"bimetallic") == 0 )
   {

     for(i=0;i<Nswap;i++)
     {
        a=M+(double)rand() / ((double)RAND_MAX / (N-M+1)+1);
        b=M+(double)rand() / ((double)RAND_MAX / (N-M+1)+1);

        if(strcmp(atom[a].Symbol.c_str(),atom[b].Symbol.c_str()) != 0 )
        {
           aux=atom[a].Symbol;
           aux2=atom[b].Symbol;

           atom[b].Symbol=aux;
           atom[a].Symbol=aux2;

        }
        else
        {
        i--;
        }
     }



   }

}


/******************************   +   *******************************/
/**************** molecule3 = molecule1 + molecule2 *****************/
/********************************************************************/

Atomic_Structure operator +(Atomic_Structure Mol1, Atomic_Structure Mol2)
{
   int j;
   Atomic_Structure Joined;
   Joined.Nat=Mol1.Nat+Mol2.Nat;

   Joined.atom=new Atom[Joined.Nat+1];

   for(i=0;i<Mol1.Nat;i++)
   {
      Joined.atom[i].Symbol=Mol1.atom[i].Symbol;
      for(j=0;j<3;j++)
      {
         Joined.atom[i].x[j]=Mol1.atom[i].x[j];
      }
   }

   for(i=0;i<Mol2.Nat;i++)
   {
   Joined.atom[i+Mol1.Nat].Symbol=Mol2.atom[i].Symbol;
      for(j=0;j<3;j++)
      {
         Joined.atom[i+Mol1.Nat].x[j]=Mol2.atom[i].x[j];
      }
   }
   return Joined;
}




/********************* Function VASP_to_xyz  ************************/
/************** vasp_to_xyz(inputfile,outputfilet) ******************/
/***************** vasp_to_xyz(argv[1],argv[2]) *********************/

void VASP_to_xyz(string inputfile, string outputfile)
{

  struct Atom
  {
     string Symbol;
     float x[3];
  };

   struct Atom xyz[500];
   struct Atom poscar[500];

   int i, j, l, m;
   string Symbol[500];
   int N_Symbol[500];
   float Factor;
   float M[3][3];
   float Mi[3][3];
   float suma;
   int Nat;
   float s;
   int Ntyp;
   int Cartesian;
   int Sel;
   int Direct;
/*************** Selectivedynamics **********************/
   string command;
          command="cat "+inputfile+"  >> aux";
   system(command.c_str());
   system("grep \"elective\" aux | wc -l >>sel ; rm aux ");

   ifstream po("sel");
            po>>Sel;
            po.close();
   system("rm sel");
   command.clear();

   if(Sel==1)//Si hay selective dynamics)
   {
      command="tail -$(($(grep -A 500 \"elective\" "+inputfile+" | wc -l )-2)) ";
      command+= inputfile+" | grep .  >> selective  ";

      ofstream tg("commandi");
               tg<<command<<endl;
               tg.close();
      system("chmod +x commandi ");
      command.clear();

      system(" ./commandi ");
      system("rm commandi ");
      system("Nat=$(cat selective | wc -l  ) " );

      system(" awk '{print $4 \" \" $5 \" \" $6}' selective | tr 'T' '1' |tr   'F' '0' >>selectivedynamicsaux ");
      system( " echo \" \" >> selectivedynamics ; echo \" \" >>selectivedynamics ; cat selectivedynamicsaux >> selectivedynamics");
      system("rm selective selectivedynamicsaux");
   }
   else{}  //Creo que esto no importa
/************************************************/
   command = "grep \"irect\" "+inputfile+" | wc -l >> Direct";
   system(command.c_str());

   ifstream x("Direct");
            x>>Direct;
            x.close();
   system("rm Direct");
   command.clear();

   command= "grep \"artesian\" "+inputfile+" | wc -l >> Cartesian";
   system(command.c_str());
   command.clear();
   ifstream xi("Cartesian");
            xi>>Cartesian;
            xi.close();
   system("rm Cartesian");

   if(Direct==1)
   {
      command="head -5 "+ inputfile+ " | tail -3 >> Matriz";
      system(command.c_str());
      command.clear();
      command="head -2 "+inputfile+ " | tail -1 >> Factor";
      system(command.c_str());
      command.clear();
      command= "tail -$(($(grep -A 500 \"irect\" "+inputfile+" | wc -l )-1)) ";
      command+= inputfile+ " | awk '{ print $1 \"  \" $2 \"  \" $3  }' | grep . >> Posiciones ";
      system(command.c_str());
      command.clear();
      command= "head -6 "+inputfile+"  | tail -1 >> aux";
      system(command.c_str());
      command.clear();
      command= "head -7 "+inputfile+"  | tail -1 >> aux2";
      system(command.c_str());
      command.clear();
      system("tr -s '[:blank:]' '\n' < aux >> Simbolos1");
      system("tr -s '[:blank:]' '\n' < aux2 >> Numeros");
      system("grep  . Simbolos1 >> Simbolos");
      system("rm Simbolos1");
      system("cat Simbolos | wc -l >>  Ntyp");

      ifstream f("Matriz");
      ifstream g("Factor");
      ifstream h("Posiciones");
      ifstream q("Simbolos");
      ifstream r("Numeros");
      ifstream p("Ntyp");
      ofstream o(outputfile);

      p>>Ntyp;
      suma=0;

      for(i=0;i<Ntyp;i++)  //Aca imagino podriamos hacer el feof c++
      {
         q>>Symbol[i];
         r>>N_Symbol[i];
         suma=suma+N_Symbol[i];
      }
      Nat=suma;
      l=0;
      for(i=0;i<Ntyp;i++) //Reconoce que tipo de atomo es cada vector posicion
      {
         for(j=0;j<N_Symbol[i];j++)
         {
            xyz[l].Symbol= Symbol[i];
            l=l+1;
         }
      }

      for(i=0;i<3;i++) //Lee los elementos de matriz
      {
         f>>Mi[i][0]>> Mi[i][1]>>Mi[i][2];
      }

      g>>Factor;  //Lee el factor de escala

      for(i=0;i<3;i++)//Multiplica el factor de escala a la matriz
      {
         for(j=0;j<3;j++)
         {
            M[i][j]=Factor*M[j][i];
         }
      }
      for(i=0;i<Nat;i++) //Extrae las coordenadas de Positions
      {
         h>>poscar[i].x[0]>>poscar[i].x[1]>>poscar[i].x[2];
      }

      for(i=0;i<Nat;i++)//Aplica la matriz  de cambio a cada atomo
      {
         for(l=0;l<3;l++)
         {
            for(m=0;m<3;m++)
               {
                  M[l][m]=poscar[i].x[l]*Mi[l][m];
               }
         }
      for(m=0;m<3;m++)
      {
         s=0;
         for(l=0;l<3;l++)
         {
            s=s+M[l][m];
         }
         xyz[i].x[m]=s;
      }
   }
   system("rm Matriz Factor Posiciones Ntyp Simbolos Numeros aux aux2"); //Borra archivos residuales
   o<<Nat<<endl<<endl;

   for(i=0;i<Nat;i++) //Imprime las lineas de atomos
   {
      o<<xyz[i].Symbol<<" "<< xyz[i].x[0]<<" "<< xyz[i].x[1]<<" "<< xyz[i].x[2]<<endl;
   }
   o.close();

   }//Cierra if/********************************************************/
   else
   {
      if(Cartesian==1)
      {
         command.clear();
         command="tail -$(($(grep -A 500 \"artesian\" "+inputfile+" | wc -l )-1)) "+inputfile+" | awk '{ print $1 \"  \" $2 \"  \" $3  }' | grep . >> Posiciones ";
         system(command.c_str());
         command.clear();
         command= "head -6 "+ inputfile+ "  | tail -1 >> aux";
         system(command.c_str());
         command.clear();
         command= "head -7 "+ inputfile+"  | tail -1 >> aux2";
         system(command.c_str());
         command.clear();
         command="tr -s '[:blank:]' '\n' < aux >> Simbolos1";
         system(command.c_str());
         command.clear();
         command="tr -s '[:blank:]' '\n' < aux2 >> Numeros";
         system(command.c_str());
         command.clear();
         system("grep  . Simbolos1 >> Simbolos");
         system("rm Simbolos1");
         system("cat Simbolos | wc -l >>  Ntyp");//Guarda el numero de especies como archivo


         ifstream hi("Posiciones");
         ifstream qi("Simbolos");
         ifstream ri("Numeros");
         ifstream pi("Ntyp");
         ofstream oi(outputfile);

         pi>>Ntyp;

         for(i=0;i<Ntyp;i++)
         {
            qi>>Symbol[i];
            ri>>N_Symbol[i];
            suma=suma+N_Symbol[i];
         }
         Nat=suma;
         l=0;
         for(i=0;i<Ntyp;i++) //Reconoce que tipo de atomo es cada vector posicion
         {
            for(j=0;j<N_Symbol[i];j++)
            {
               xyz[l].Symbol=Symbol[i];
               l=l+1;
            }
         }
         for(i=0;i<Nat;i++) //Extrae las coordenadas de Positions
         {
            hi>>poscar[i].x[0]>>poscar[i].x[1]>>poscar[i].x[2];
         }
         system("rm  Posiciones Ntyp Simbolos Numeros aux aux2"); //Borra archivos residuales
         oi<<Nat<<endl<<endl;
         for(i=0;i<Nat;i++) //Imprime las lineas de atomos
         {
            oi<<xyz[i].Symbol<<" "<< poscar[i].x[0]<<" "<< poscar[i].x[1]<<" "<< poscar[i].x[2]<<endl;
         }
         oi.clear();
      }
   }

   if (Sel==1)
   {
      command="cat "+outputfile+" >> aux ;  rm "+outputfile+" ; paste aux selectivedynamics  >> "+outputfile;
      system(command.c_str());
      system("rm  aux selectivedynamics");
   }
}//Cierra funcion


/********************* Function xyz_to_VASP  ************************/
/** xyz_to_vasp(InFile,OutFile,Title,ScaleFactor,Matrix(3x3Array) ***/
/******float M[3][3]; xyz_to_vasp(mol.xyz,POSCAR,"Title",1,M);*******/

void xyz_to_VASP(string inputfile,string outputfile,string Titulo, float Factor, float M[3][3])  //Checa la matriz
{
   int i,j;
   int Selective;
   string command;
   ofstream hj(outputfile);
   hj<<Titulo<<endl;
   hj<<Factor<<endl;

   for(i=0;i<3;i++)           //Imprime Matriz de vectores
   {
      for(j=0;j<3;j++)
      {
         hj<<M[i][j]<<"  ";
      }
      hj<<endl;
   }
   hj.close();

   command="cat "+inputfile+" >> selective.txt";
   system(command.c_str());
//   system("head -$(($(head -1 selective.txt )+2)) selective.txt | tail -$(head -1 selective.txt) | awk '{print $5}'");
   command.clear();
   command="awk '{print $5 }'  selective.txt  | grep . | wc -l >> sel";
   system(command.c_str());

   ifstream iy("sel");
   iy>>Selective;
   iy.close();
   system("rm sel");
   if(Selective>0)
   {
      system("head -$(($(head -1 selective.txt )+2)) selective.txt | tail -$(head -1 selective.txt) | awk '{print $5  \" \" $6 \" \" $7 }'  | tr '0' 'T'  >> aux");
      system("rm selective.txt ; cat aux >> selective.txt ; rm aux ");
   }
/////////////////////////////////||||||||||||

   ofstream pl("command");
   pl<< "#!/bin/bash"<<endl;
   pl<<"N=$(head -1 $1)"<<endl;
   pl<<"head -$(($N+2)) $1 | tail -$N |  awk '{print $1 }' | sort -u >> aux"<<endl;
   pl<<"M=$(wc -l aux | awk '{ print $1 }')"<<endl;
   pl<<"for ((i=1;i<$(($M+1));i++))"<<endl;
   pl<<"do"<<endl;
   pl<<"echo -n \"$(head -$i aux | tail -1  )  \" >>$2"<<endl;
   pl<<"done"<<endl;
   pl<<"echo \" \" >>$2"<<endl;

   pl<<"for ((i=1;i<$(($M+1));i++))"<<endl;
   pl<<"do"<<endl;
   pl<<"echo -n \"$( grep -w \"$(head -$i aux | tail -1 )\" $1 | wc -l )  \">> $2"<<endl;
   pl<<"done"<<endl;
   pl<<"echo \" \" >>$2"<<endl;

   if(Selective>0) //Lo que significa que hay al menos una linea seleccionada
   {
      pl<<"echo \"Selective dynamics\" >>$2 "<<endl;
   }
   pl<<"echo \"Cartesian\">>$2 "<<endl;


   pl<<"rm aux"<<endl;
   pl<<"head -$(($N+2)) $1 | tail -$N |  sort -u | awk '{print $2 \"  \"$3 \"  \" $4  }'  >> aux1"<<endl;
   pl<<"head -$(($N+2)) $1 | tail -$N |  sort -u | awk '{print $5 \"  \"$6 \"  \" $7  }'  >> aux2"<<endl;

   if(Selective>0) //Lo que significa que hay al menos una linea seleccionada
   {
      pl<<"Nl=$(cat aux2 | wc -l ) \n for((i=1;i<$(($N+1)); i++)) \n do \n cont=$(head -$i aux2 | tail -1 | tr '0' 'F' | tr '1' 'T' ) "<<endl;
      pl<<" ki=$(echo $cont | wc -c) \n if [ $ki -gt 1 ] \n then echo \"$cont\" >>aux3 \n else \n echo \" T   T   T \" >> aux3  "<<endl;
      pl<<"fi \n done \n paste aux1 aux3 >> $2 \n rm aux1 aux2   aux3 "<<endl;
   }
   else
   {
      pl<<"paste aux1 aux2 >>$2 ; rm aux2 aux1  "<<endl;
   }
   pl.close();
   system("chmod +x command");

   command.clear();
   command="./command  "+inputfile+ " "+outputfile;
   system(command.c_str());
   system("rm selective.txt command");
}


/**************************** read_VASP *****************************/
/********************* Molecule  molecule_name; *********************/
/**************** molecule_name.read_VASP(argv[1]); *****************/
/********************************************************************/

void Atomic_Structure::read_VASP(string file)
{
   VASP_to_xyz(file, "salida");
      float x,y,z;
      string Symbol;
      string command;
      system("head -1 salida >> Nat");
      ifstream Nat_file("Nat");
               Nat_file >> Nat;
               Nat_file.close();
               cout<<Nat<<endl;

      system("rm Nat");
             command.clear();
             command="head -";
             command+=to_string(Nat+2);
             command+=" salida | tail -";
             command+=to_string(Nat);
             command+=" >> coordinatesAux";
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
      system("rm coordinatesAux ");

   system("rm salida");
   map<string, double> Radios;
   Radios=radii_dictionary();
   for(i=0;i<Nat;i++)
   {
      atom[i].R=assign_radii(Radios,atom[i].Symbol);
   }
}


/***************************** centroid *****************************/
/********************* molecule_name.centroid() *********************/
/************************ Cysteine.centroid *************************/
/********************************************************************/

void Molecule::centroid()
{
   float Dx, Dy, Dz, sumax=0,sumay=0,sumaz=0;
   for(i=0;i<Nat;i++)
   {
      sumax=atom[i].x[0]+sumax;
      sumay=atom[i].x[1]+sumay;
      sumaz=atom[i].x[2]+sumaz;
   }
   Dx=sumax/Nat;
   Dy=sumay/Nat;
   Dz=sumaz/Nat;

   for(i=0;i<Nat;i++)
   {
       atom[i].x[0]=atom[i].x[0]-Dx;
       atom[i].x[1]=atom[i].x[1]-Dy;
       atom[i].x[2]=atom[i].x[2]-Dz;
   }
}
/***************************** centroid *****************************/
/********************* molecule_name.centroid() *********************/
/************************ Cysteine.centroid *************************/
/********************************************************************/

void Cluster::centroid()
{
   float Dx, Dy, Dz, sumax=0,sumay=0,sumaz=0;
   for(i=0;i<Nat;i++)
   {
      sumax=atom[i].x[0]+sumax;
      sumay=atom[i].x[1]+sumay;
      sumaz=atom[i].x[2]+sumaz;
   }
   Dx=sumax/Nat;
   Dy=sumay/Nat;
   Dz=sumaz/Nat;

   for(i=0;i<Nat;i++)
   {
       atom[i].x[0]=atom[i].x[0]-Dx;
       atom[i].x[1]=atom[i].x[1]-Dy;
       atom[i].x[2]=atom[i].x[2]-Dz;
   }
}

/**************************** print_VASP ****************************/
/********************* Molecule  molecule_name; *********************/
/**************** molecule_name.print_VASP(argv[1]); ****************/
/********************************************************************/

void Atomic_Structure::print_VASP(string outputfile,string Titulo, float Factor, float M[3][3])
{
   ofstream output_file("output");
   output_file<<Nat<<endl;
   output_file<<" "<<endl;
   for(i=0;i<Nat;i++)
   {
      output_file<<atom[i].Symbol<<" "<<atom[i].x[0]<<" "<<atom[i].x[1]<<" "<<atom[i].x[2]<<endl;
   }
   output_file.close();
system("cat output");
   xyz_to_VASP("output", outputfile, Titulo, Factor, M);
   system("rm output");
}


/*

int main()
{
float Matrix[3][3]={{25,0,0},{0,25,0},{0,0,25}};
cout<<"Variables definidas: Matriz"<<endl;

Atomic_Structure Conf;
cout<<"Variables definidas conf"<<endl;

Molecule Cysteine;
cout<<"Variables definidas molecule"<<endl;

Cluster Au20;
cout<<"Variables definidas"<<endl;
Au20.rand_generator("U",20);
cout<<"cluster generado"<<endl;
Cysteine.read_VASP("POSCAR");

cout<<"moleculaomovida"<<endl;
Cysteine.rotate_Rad(3.14,0.0);
Au20.rotate_Deg(90,45);
Conf=Cysteine+Au20;
cout<<"configuracion creada"<<endl;


cout<<"==========="<<endl;
Cysteine.print_VASP("Cyst.vasp","TITULO",1.0,Matrix);
Au20.print_VASP("Clus.vasp","TITULO",1.0,Matrix);
Conf.print_VASP("CONTCARPru","TITULO",1.0,Matrix);
return 0;
}*/
