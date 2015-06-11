// Modulation: Gleeson & Axford (1968) force-field approximation /imos 2002/7/30
// Ekin=MeV, cr_density = CR density (E^2*flux) , phi=MV, key=0,1
// AWS20020730 float *cr_density (original had double *)

/*
AWS20091013
Documentation of method.

 Gleeson & Axford,  ApJ 154, 1011 (1968)
The formulae here are based on equation (16) of that paper.

 J(E)           J(E + |Z| e phi)
________   =    _________________________ 
E^2-Eo^2         (E + |Z| e phi)^2 - Eo^2        

where E = total energy of particle, Eo = rest mass energy of particle.
Rewriting this using E=ekin+m where ekin = kinetic energy of particle, m=rest mass of particle
t = ekin + |Z| e phi
 E^2              - Eo^2 =  ekin (ekin + 2m)
(E + |Z| e phi)^2 - Eo^2  =  t   (t    + 2m)

                        ekin (ekin + 2m) 
J(ekin)   =     J(t) *  ________________
                         t   (t    + 2m)

   In GALPROP we work always in Ekin per nucleon, not total Ekin of particle. Hence we have to rewrite the formulae.

   For nuclei with m=A Mp, and writing Ekin = ekin/A   to correspond to the notation in the routine

   J(Ekin) dEkin = J(ekin) dekin
   J(Ekin)       = J(ekin) dekin/dEkin = A J(ekin)

   define T also per nucleon:
   T = Ekin + |Z| e phi / A        i.e. the potential energy loss is divided between the nucleons, hence the division by A, and phi-> phi/A
   t = E    + |Z| e phi = A Ekin +|Z| e phi/A * A = AT

   J(T)         = J(t)   dt/dT         = A J(t) 



                               A Ekin (A Ekin + 2A Mp)               Ekin (Ekin + 2 Mp) 
     J(Ekin)   =     J(T) *    _______________________      = J(T)   __________________
                               A  T   (A T    + 2A Mp)                T   (T    + 2 Mp)


   Finally, since the input and output J are expressed as Ekin^2 J(Ekin) (a GALPROP convention), and the calculations are all done with this quantity,
   the final result has to be  multiplied by (Ekin/T)^2:



   Explitly:
   Ekin^2 J(Ekin) =  Ekin^2 J(T) * [Ekin (Ekin + 2Mp)] /[  T   (T    + 2Mp)]
                  =     T^2 J(T) * [Ekin (Ekin + 2Mp)] /[  T   (T    + 2Mp)]  [Ekin^2/ T^2]
         

   

*/

using namespace std; //AWS20050919
#include <cstdlib>   //AWS20050919
#include <cmath>     //AWS20050919
#include <fstream>  //AWS20050919
#include <iostream> //AWS20091009

void modulate(double *Ekin, float *cr_density, int np, int z, int a, double phi, int key)
{
   int j;
   double *density,B,C,T,y,Mp=939.; // Mp = (mp+mn)/2 ~939 MeV
   double                  Me=0.511; //AWS20091009
   density=new double[np];

   int debug = 0; //AWS20091009
 
   for(int i=0; i<np; i++)
   {
     if(a> 0) T = Ekin[i]+abs(z)*phi/a;       //AWS20091012
     if(a==0) T = Ekin[i]+abs(z)*phi        ; //AWS20091012 electrons, positrons

      for(j=0; j<np; j++) if(T<Ekin[j]) break;
      if(j==np) { density[i] = cr_density[i]; break; }
      switch(key)
      {
         case 0: // linear interpolation in energy
	    y = cr_density[j-1] +(T-Ekin[j-1]) *(cr_density[j]-cr_density[j-1]) /(Ekin[j]-Ekin[j-1]);
            break;
         case 1: // power-law interpolation
	    if(cr_density[j-1]*cr_density[j]>0.)
	    {
	       B = log(cr_density[j-1]/cr_density[j])/log(Ekin[j-1]/Ekin[j]);
               C = cr_density[j]*pow(Ekin[j],-B);
               y = C*pow(T,B);
            }
            else  // linear interpolation in energy
	       y = cr_density[j-1] +(T-Ekin[j-1]) *(cr_density[j]-cr_density[j-1]) /(Ekin[j]-Ekin[j-1]);
            break;
         default:
	    y = cr_density[j-1] +(T-Ekin[j-1]) *(cr_density[j]-cr_density[j-1]) /(Ekin[j]-Ekin[j-1]);
       }                                                    // remove last factor if working with flux 
    
      density[i]= y * Ekin[i] *(Ekin[i]+2*Mp) /T /(T+2*Mp)    *pow(Ekin[i]/T,2);

      if(a==0)                                                                   // electrons, positrons AWS20091009
      density[i]= y * Ekin[i] *(Ekin[i]+2*Me) /T /(T+2*Me)    *pow(Ekin[i]/T,2); // electrons, positrons AWS20091009

      if(debug==1&&a==0)
	{
	  cout<<"modulate: z="<<z<<" a="<<a<<" Ekin="<<Ekin[i]<<" T="<<T<<"  cr_density="<<cr_density[i]<<" density="<<density[i]<<endl;
	}

   }
   for(j=0;j<np;j++) cr_density[j] = density[j];

   delete[] density;
}

/*
#include <fstream.h>
void modulate(double*,float*,int,int,int,double,int);
main()
{
   int i=0,z=6,a=12,key;
   double Ekin[100],phi=500.,tmp[100],tmp1;
   float cr_density[100],cr[100];
   ifstream data;
   data.open("42.2_111331");

   while(!data.eof())
   {
      data>>Ekin[i]>>cr_density[i]>>tmp[i]>>tmp1>>tmp1;
      cout<<"input> "<<Ekin[i]<<" "<<cr_density[i]<<" "<<tmp[i]<<endl;
      i++;
   }
   cout<<endl<<"output:"<<endl;
   for(int j=0;j<i-1;j++) cr[j] = cr_density[j];
   modulate(Ekin, cr_density, i-1, z, a, phi, 0);
   modulate(Ekin, cr        , i-1, z, a, phi, 1);
  for(int j=0;j<i-1;j++)
   {
      cout<<Ekin[j]<<" "<<cr_density[j]<<" "<<cr[j]<<" "<<tmp[j]<<endl;
   }
}

*/
