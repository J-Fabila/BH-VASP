
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_ RandomGenerator.cpp _/_/_/_/_/_/_/_/_/_/_/_/
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

using namespace std;

/********* * * * * * *  *  *   *   *   *  *  *  * * * ***************/

int Nat;
int i,j,k;

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
      void read_VASP(string);
      void centroid();
      double x_min(); double x_max();
      double y_min(); double y_max();
      double z_min(); double z_max();
      bool fit_in(float, float, float, float, float, float);

};

using Crystal=Molecule;
using Cluster=Molecule;
using Atomic_Structure=Molecule;


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

/**************************** print_xyz *****************************/
/********************* Molecule  molecule_name; *********************/
/**************** molecule_name.print_xyz(argv[1]); *****************/
/********************************************************************/

void Molecule::print_xyz(string outputfile)
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

/*
         atom[i].x[0]=Mx + ( rand() % ( Mx - Nx + 1 ) );
         atom[i].x[1]=My + ( rand() % ( My - Ny + 1 ) );
         atom[i].x[2]=Mz + ( rand() % ( Mz - Nz + 1 ) );
*/
         atom[i].x[0]=Nx+ ((double)rand())/((double)RAND_MAX )* (Mx-Nx);
         atom[i].x[1]=Ny+ ((double)rand())/((double)RAND_MAX )* (My-Ny);
         atom[i].x[2]=Nz+ ((double)rand())/((double)RAND_MAX )* (Mz-Nz);

/*       atom[i].x[1]=Ny+(double)rand()/((double)RAND_MAX/(Ny-My+1)+1);
         atom[i].x[2]=Nz+(double)rand()/((double)RAND_MAX/(Nz-Mz+1)+1);

         atom[i].x[0]=Nx+(double)rand()/((double)RAND_MAX/(Nx-Mx+1)+1);
         atom[i].x[1]=Ny+(double)rand()/((double)RAND_MAX/(Ny-My+1)+1);
         atom[i].x[2]=Nz+(double)rand()/((double)RAND_MAX/(Nz-Mz+1)+1);
*/
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


/**************************** read_VASP *****************************/
/********************* Molecule  molecule_name; *********************/
/**************** molecule_name.read_VASP(argv[1]); *****************/
/********************************************************************/

void Molecule::read_VASP(string file)
{
   VASP_to_xyz(file, "salida");
      float x,y,z;
      string Symbol;
      string command;
      system("head -1 salida >> Nat");
      ifstream Nat_file("Nat");
               Nat_file >> Nat;
               Nat_file.close();

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

float random_number(float min,float max)
{
   float random = min + ((float)rand())/((float)RAND_MAX ) * (max-min);
   return random;
}

int main(int argc, char *argv[])
{
double _x, _y, _z;
Cluster clus;
if(argc==11)
{
  clus.random_generator(argv[1],stoi(argv[2]));
  bool acep=clus.fit_in(stoi(argv[3]),stoi(argv[4]),stoi(argv[5]),stoi(argv[6]),stoi(argv[9]),stoi(argv[10]));
  while (acep==false)
  {
     clus.~Molecule();
     clus.random_generator(argv[1],stoi(argv[2]));
     acep=clus.fit_in(stoi(argv[3]),stoi(argv[4]),stoi(argv[5]),stoi(argv[6]),stoi(argv[9]),stoi(argv[10]));
   }
  _x= stoi(argv[3]) + ((double)rand())/((double)RAND_MAX )* (stoi(argv[4])-stoi(argv[3]));
  _y= stoi(argv[5]) + ((double)rand())/((double)RAND_MAX )* (stoi(argv[6])-stoi(argv[5]));
  _z= stoi(argv[9]);

}
else
{
if(argc==13)
{

clus.random_generator(argv[1],stoi(argv[2]),argv[3],stoi(argv[4]));
bool acep=clus.fit_in(stoi(argv[5]),stoi(argv[6]),stoi(argv[7]),stoi(argv[8]),stoi(argv[11]),stoi(argv[12]));
while (acep==false)
{
   clus.~Molecule();
   clus.random_generator(argv[1],stoi(argv[2]),argv[3],stoi(argv[4]));
   acep=clus.fit_in(stoi(argv[5]),stoi(argv[6]),stoi(argv[7]),stoi(argv[8]),stoi(argv[11]),stoi(argv[12]));
 }
_x= stoi(argv[5]) + ((double)rand())/((double)RAND_MAX )* (stoi(argv[6])-stoi(argv[5]));
_y= stoi(argv[7]) + ((double)rand())/((double)RAND_MAX )* (stoi(argv[8])-stoi(argv[7]));
_z= stoi(argv[11]);
}
}
double zmin=clus.z_min();
clus.move( _x , _y , -zmin + _z );

  clus.print_xyz("ClusterGenerated.xyz");
  float rnd=random_number(10,20);
  return 0;
}
