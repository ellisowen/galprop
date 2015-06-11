

// g++ Bremsstrahlung.cc inter_cc.cc sim.cc test_Bremsstrahlung.cc ../ProductionMatrices/CR_Spectrum.cc -I../ProductionMatrices -o test_Bremsstrahlung

//  rewritten  from original fortran in GALPROP public version (see documentation appended to this file).
//  tested against original version, mainly agrees but some cases still do not.

// A. Strong, 16 Apr 2012
//            20 Jun 2012 : pow(int,int) -> pow(int,float) for 3 cases to satisfy some compilers. Corrected compile instructions: inter_cc.cc 

////////////////////////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace FPHI{double  delta; int IZ; int Ne;} // use namespace to share with  FPHI1  FPHI2
using namespace  FPHI;

double   FPHI1(double Q), FPHI2(double Q);
 
double sim(double, double, double, double, double, double(*)(double));

void inter(double*, double*, int, double, double&);


#include<iostream>
#include<cmath>
#include"Bremsstrahlung.h"
#include"CR_Spectrum.h"

/////////////////////////////////////////////////////////////

void  Bremsstrahlung::  init(string name_, valarray<double> parameters_, int debug_)
{
 cout<<">> Bremsstrahlung::init "<<endl ;

 name=name_;
 parameters=parameters_;
 debug=debug_;
 initialized=1;

 cout<<"<< Bremsstrahlung::init"<<endl ;
 return;

};

/////////////////////////////////////////////////////////////

double Bremsstrahlung::dSdE(double Ee, double Egam, int IZ1, int Ne1 )
{
  // dSdE = differential cross-secion in barn/MeV
  // Ee   =  electron energy in MeV (NB total energy)
  // Egam = gamma-ray energy in MeV
  // IZ1  = nucleus charge
  // Ne1  = 0 for unscreened target nucleus, 1 for 1-electron atom, 2 for 2-electron atom
  // For details see original documentation appended to this file

  if(debug==1) cout<<">> Bremsstrahlung::dSdE" ;

  
  if(debug==1) cout<< "  input parameters  Ee="<<Ee<<" MeV Egam="<<Egam<< " MeV"
		   <<" charge=IZ1="<<IZ1<<" number of electrons= Ne1="<<Ne1<<endl;

  // He formfactors  size=12 to allow indices 1-11 as in original fortran version
  double DD [12]={-999,   0.,     0.01,   0.02,   0.05,   0.1,   0.2,    0.5,    1.,     2.,     5.,    10.  };
  double PH1[12]={-999, 134.60, 133.85, 133.11, 130.86, 127.17,120.35, 104.60,  89.94,  74.19,  54.26,  40.94};
  double PH2[12]={-999, 131.40, 130.51, 130.33, 129.26, 126.76,120.80, 105.21,  89.46,  73.03,  51.84,  37.24}; 
  


  double DD1[12]; // size=12 to allow indices as in original fortran version
  double phi1,phi2,p1,p2;

  double dSdK = 0.;
  double mc2 = 0.511; //    MeV, electron rest mass

  double E0 = Ee;     //   original notation, but argument Ee for compatibility with other tools
  double E1 = E0-Egam;//   total energy of  electron after scattering

  


  if(E0 <= 1.02*mc2 || E1 <= mc2)
  {
   if(debug==1)cout<<"Bremsstrahlung dSdE part 0: E0 <= 1.02*mc2 || E1 <= mc2  "
		   <<   "  Ee="<<Ee<<" MeV Egam="<<Egam<< " MeV"<<" E1="<<E1<<" T0=E0-mc2="<<E0-mc2
		   <<" IZ1="<<IZ1<<" Ne1="<<Ne1
		   <<" dSdK="<<dSdK<<endl;
   return dSdK;
  }

         

 IZ    = IZ1;                 // in namespace, see top
 Ne    = Ne1;                 // in namespace, see top
 delta = Egam*mc2/(2.*E0*E1); // in namespace, see top

 double Pi = 3.14159265359;
 double Af = 7.2974e-3;          //       fine structure constant
 double R02 =pow( 2.81794e-1,2); //  classical electron radius square (barn)
 double EB1 = 0.07;              //  MeV,kinetic energy; the upper bound for NonRelat.approx.
 double EB2 = 2.;                //  MeV, -"-; the lower bound for HE approx. (nonscreened)



 double Ef = 0.;
 double FS = 0.;
 double FS90 = 0.90;
 double T0 = E0-mc2 ;                     //   Ekin
 double P0 = mc2* sqrt(pow(E0/mc2,2)-1.); //  initial electron momentum
 double P1 = mc2* sqrt(pow(E1/mc2,2)-1.); //    final electron momentum
 double beta0 = P0/E0;
 double beta1 = P1/E1;

  // Fano-Sauter high-frequency limit (Egam -> E0-mc2, KM59,p.933) 

      if(Egam/T0 >=   FS90)
	FS =4.*Pi*pow(IZ,3.)*R02*pow(Af,2)/T0*beta0*E0/mc2/pow(E0/mc2-1.,2) 
            *(4./3. +E0/mc2*(E0/mc2-2.  )/(E0/mc2+1.  )
	    *(1.- log(1.  +2.  *beta0/(1.  -beta0))*pow(mc2/E0,2)/2./beta0));

  // Elwert factor (Ef) is valid if   Z*Af*(1./beta1-1./beta0)<<1 

      if(T0 <    EB2 &&    IZ*Af*(1./beta1-1./beta0) <=   0.5)
         Ef = beta0/beta1 *(1.  - exp(-2.  *Pi*Af*IZ/beta0))
	   /(1.  - exp(-2.  *Pi*Af*IZ/beta1));


   if(debug==1)cout<<"Bremsstrahlung dSdE before part 1: E0 <= 1.02*mc2 || E1 <= mc2  "
		   << "  Ee="<<Ee<<" MeV Egam="<<Egam<< " MeV"<<" E1="<<E1<<" T0=E0-mc2="<<E0-mc2<<" EB1="<<EB1
		   <<" IZ1="<<IZ1<<" Ne1="<<Ne1
                   <<" FS="<<FS<<" Ef="<<Ef
		   <<endl;

  // NONRELATIVISTIC ENERGIES ( T0 < 0.07 MeV) 
  // Born approx.(KM59,p.925,3BNa) with Elwert factor (p.931)
  // Born approximation: 2*Pi*Z*Af/beta{0,1}<<1; 
  // nonscreened: 2*delta/[Af*Z^(1/3)] >> 1;

      if(T0 <    EB1) 
      {
	  int A = IZ*IZ ;                        // no electron contribution
          dSdK  = Ef*R02*Af*A*16./3. /Egam 
	         *pow(mc2/P0,2)*log(1.+2.*P1/(P0-P1));

	  if(dSdK < FS) dSdK = (FS+dSdK)/2.;     // 1/2 Fano-Sauter limit


          if(debug==1)cout<<"Bremsstrahlung dSdE part 1: T0 <    EB1: "
			  <<"  Ee="<<Ee<<" MeV Egam="<<Egam<< " MeV"<<" E1="<<E1<<" T0=E0-mc2="<<E0-mc2<<" EB1="<<EB1
		          <<" IZ1="<<IZ1<<" Ne1="<<Ne1
			  <<" dSdK="<<dSdK<<endl;
	  return dSdK;
       }

   // LOW ENERGY (T0 < 2. MeV) OR delta/(2.*Af*Z) > 4. 
   // Born approximation (KM59,p.925,3BN): 2*Pi*Z*Af/beta{0,1}<<1;
   // nonscreened: 2*delta/Af >> 1

      if(T0 <    EB2 ||   delta/(2.*Af*IZ) >=   4.)     
      {
	

	 double eps0 =  log(1.  +2.  *P0/(E0-P0));// =log((E0+P0)/(E0-P0))
	 double eps1 =  log(1.  +2.  *P1/(E1-P1));// =log((E1+P1)/(E1-P1))
	 double AL   = 2.  * log((E0*E1+P0*P1-pow(mc2,2))/(Egam*mc2));
	 double T1   = 2.  *E0*E1*(pow(P0/P1,2)+1.  )/pow(P0,2);
	 double T2   = pow(mc2,2)*(eps0*E1/pow(P0,3)+eps1*E0/pow(P1,3)-eps0*eps1/P0/P1);
	 double T10  = 8.  /3.*(E0*E1)/(P0*P1);
	 double T11  = pow(Egam,2)*(pow(E0/P0*E1/P1,2)+1.  )/(P0*P1);

	 double T12  = pow(mc2,2) *Egam/(2.  *P0*P1)
                      *(
                        eps0*(E0/P0*E1/P0+1.  )/P0 -eps1*(E0/P1*E1/P1+1.  )/P1
                        +2.  *Egam*E0/P0*E1/P1/(P0*P1) );

         int A = IZ*IZ +Ne ;                                          //  correction for atomic electrons

         if(T0 <   EB2)
	   A = ( IZ*IZ +Ne*(1.  - exp(-(T0-EB1)/9./EB1)) )            //  screening
	      *(1.-0.3* exp(-Egam/0.033)) *Ef;                        //  left & right correction

         dSdK = R02*Af*A/Egam*P1/P0*(4.  /3.-T1+T2+AL*(T10+T11+T12));

         if(dSdK <    FS) dSdK = (FS+dSdK)/2.;                        //  1/2 Fano-Sauter limit


          if(debug==1)cout<<"Bremsstrahlung dSdE part 2: T0 <    EB2 ||   delta/(2.*Af*IZ) >=   4. : "
			  <<   "  Ee="<<Ee<<" MeV Egam="<<Egam<< " MeV"
		          <<" IZ1="<<IZ1<<" Ne1="<<Ne1
			  <<" dSdK="<<dSdK<<endl;
	 return dSdK;                   
	}
      



  // HIGH  ENERGY (T0 > 2. MeV) AND delta/(2.*Af*Z) < 4. 

  // Schiff approx. for neutral atoms heavier then He (KM59,p.925,3BSe) 

  if(IZ > 2 && IZ ==   Ne)      
  {
    double b = pow(IZ, 1./3.)/111.  /delta;
    double OM = 1.  /pow(delta,2)/(1.  +b*b);

         dSdK = 2.  *IZ*IZ*R02*Af/Egam *( 
		(1.  +pow(E1/E0,2)-2.  /3.*E1/E0)
               *( log(OM)+1.  -2.  /b* atan(b))
               + E1/E0*(2.  /b/b* log(1.  +b*b)
		        +4./3.*(2.  -b*b)/pow(b,3) * atan(b)-8.  /3./b/b+2.  /9.));

         if(debug==1)cout<<"Bremsstrahlung dSdE part 3: Schiff approx for neutral atoms >He: IZ > 2 && IZ ==   Ne  : "
			  <<   "  Ee="<<Ee<<" MeV Egam="<<Egam<< " MeV"
		          <<" IZ1="<<IZ1<<" Ne1="<<Ne1
			  <<" dSdK="<<dSdK<<endl;

    return dSdK;
  }

  // arbitrary screening (G69, BG70; KM59,p.925,3BSb)

  int A = IZ*IZ +Ne ;                                    // correction for atomic electrons
  double phiu = -log(delta)-1.  /2.  ;

  // Hartree-Fock approximation for neutral He atoms

      if(IZ ==   2 &&    IZ ==   Ne) 
      {
	  DD1[1] =  log(1.e-3);
	  for(int i=2;i<=11;i++)
	   DD1[i] =  log(DD[i]);
       

         if(delta/(2.  *Af) >=   DD[2]) 
	 {
            if(delta/(2.  *Af) <=   DD[11]) 
	    {
	       inter(DD1,PH1,11, log(delta/(2.  *Af)),phi1);
	       inter(DD1,PH2,11, log(delta/(2.  *Af)),phi2);
	    }

            else
	    {
		phi1 = 4.*A *phiu;//    asymptotics
		phi2 = 4.*A *phiu;
	    }

	  }

	  else
	  {
		   phi1 = PH1[1];
		   phi2 = PH2[1];
	  }

	    phi1 = phi1/4.;
	    phi2 = phi2/4.;


                if(debug==1)cout<<"Bremsstrahlung dSdE part 4: Hartree-Fock IZ ==   2 &&    IZ ==   Ne : "              
			  <<   "  Ee="<<Ee<<" MeV Egam="<<Egam<< " MeV"
		          <<" IZ1="<<IZ1<<" Ne1="<<Ne1
			  <<endl;

      } //(IZ ==   2 &&    IZ ==   Ne)

	 else
	 {
	      if(Ne ==   0)                   //   UNSHIELDED charge
	      {
		  phi1 = A*phiu;
		  phi2 = A*phiu;

                if(debug==1)cout<<"Bremsstrahlung dSdE part 5 :  Ne==0   UNSHIELDED charge : "              
			  <<"  Ee="<<Ee<<" MeV Egam="<<Egam<< " MeV"
		          <<" IZ1="<<IZ1<<" Ne1="<<Ne1
			  <<endl;
	      }
	      else
	      {
                // H-like atoms & Hylleraas-1 approximation for He-like atoms

		  if(delta/(2.  *Af) <=   1.e-4) delta = 1.e-4*(2.  *Af);

		  //		  SIM1(1.  ,delta,1.e-3,1.e-5,1.e-4,FPHI1,p1);
		  //		  SIM1(1.  ,delta,1.e-3,1.e-5,1.e-4,FPHI2,p2);


		  // NB limits inverted and then sign reversed in formula for phi1,phi2: why ?

         	  p1=	  sim (1.  ,delta,1.e-3,1.e-5,1.e-4,FPHI1   );
	          p2=	  sim (1.  ,delta,1.e-3,1.e-5,1.e-4,FPHI2   );

	          if(debug==1)cout<<" output from sim: p1="<<p1<<" p2="<<p2<<endl;

		  phi1 = 2.*IZ*(-p1       +1.+(Ne-1.)/IZ)  +pow(IZ-Ne,2.)*phiu;
		  phi2 = 2.*IZ*(-p2+5./6.*(1.+(Ne-1.)/IZ) )+pow(IZ-Ne,2.)*phiu;

                  if(debug==1)cout<<"Bremsstrahlung dSdE part 6 : H-like and He-like atoms : "              
			   <<"  Ee="<<Ee<<" MeV Egam="<<Egam<< " MeV"
		           <<" IZ1="<<IZ1<<" Ne1="<<Ne1
			   <<endl;

	      }   

    

	 }
      


	  if(phi1 >    A*phiu) phi1 = A *phiu; // asymptotics
	  if(phi2 >    A*phiu) phi2 = A *phiu;

	  dSdK =4.*R02*Af/Egam *(( 1.  +pow(E1/E0,2))*phi1 -2./3.*E1/E0*phi2);
	  if(dSdK <    FS) dSdK = (FS+dSdK)/2.;//   !# 1/2 Fano-Sauter limit


          if(debug==1)cout<<"Bremsstrahlung dSdE part 7: end of routine : "              
			  <<   "  Ee="<<Ee<<" MeV Egam="<<Egam<< " MeV"
		          <<" IZ1="<<IZ1<<" Ne1="<<Ne1
			  <<" dSdK="<<dSdK<<endl;

	  return dSdK;
      


  };




/////////////////////////////////////////////////////////////////////////////////////////////

double   FPHI1(double Q) 
{
  //c**********************************************************************
  //c used for calculation of PHI1 function (G69,p.74,76)
  //c          *** I.Moskalenko (MPE,Garching) *** version of 23.09.1997 **
  //c**********************************************************************
  //      implicit real*8 (a-h,o-z)
  // this is now in global namespace:
  //      common/delta/ delta, IZ, Ne

  int debug=0;
  if(debug==1)cout<<"FPHI1: Q= "<<Q<<endl;


  double FPHI1_;
  double AZ,FZ;
  double Af = 7.2974e-3;    //       fine structure constant

       if(Ne == 1) //      ! one-electron atoms
       {
	   AZ = 1.  /(2.  *Af*IZ);
	   FZ = 1.  /pow(1.  +pow(AZ*Q,2), 2); // form-factor for one-electron atom
	   FPHI1_ = pow(Q-delta,2)/pow(Q,3) *(1.  -FZ);
         return FPHI1_  ;
       }

  //                    two-electron atoms (Hylleraas-1 function)

      AZ = 1.  /(2.  *Af*(IZ-5.  /16.  ));
      FZ = 1.  /pow(1. +pow(AZ*Q,2),2);     //  form-factor for two-electron atom

             FPHI1_ = pow(Q-delta,2)/pow(Q,3) *( 2.*(1.-FZ) -(1.-pow(FZ,2))/IZ );
      return FPHI1_  ;
};
      

///////////////////////////////////////////////////////////////////////	

      double  FPHI2(double Q)
{ 
  //c***********************************************************************
  //c used for calculation of PHI2 function (G69,p.74,76)
  //c          *** I.Moskalenko (MPE,Garching) *** version of 23.09.1997 ***
  //c***********************************************************************
  //      implicit real*8 (a-h,o-z)
  // this is now in global namespace:
  //     common/delta/ delta, IZ, Ne

  int debug=0;
  if(debug==1)cout<<"FPHI2: Q= "<<Q<<endl;

  double  FPHI2_;
  double AZ,FZ;
  double Af = 7.2974e-3;//        fine structure constant

  if(Ne ==   1)     //       one-electron atoms
  {
	  AZ = 1.  /(2.  *Af*IZ);
	  FZ = 1.  /pow(1.  +pow(AZ*Q,2),2);//! atomic form-factor
	  FPHI2_=(pow(Q,3)+(3.e-6   * log(Q/delta))*pow(delta,2)*Q 
		  -4.  *pow(delta,3))/pow(Q,4) *(1.  -FZ);
         return  FPHI2_;
  }  
  
      //                    two-electron atoms (Hylleraas-1 function)

      AZ = 1.  /(2.  *Af*(IZ-5.  /16.  ));
      FZ = 1.  /pow(1.  +pow(AZ*Q,2),2);//! form-factor for two-electron atom

      FPHI2_=(pow(Q,3)+(3.e-6  * log(Q/delta))*pow(delta,2)*Q 
	      -4.  * pow(delta,3) )/pow(Q,4) *(2.*(1.  -FZ) -(1.  -pow(FZ,2))/IZ);

      return FPHI2_;

		      };

///////////////////////////////////////////////////////////////////////////

void  Bremsstrahlung:: bremss_emiss(valarray<double> Ee,   valarray<double> Egamma,
                                    valarray<double> fluxe,
                                    int IZ1, int Ne,
                                    valarray<double> &emiss,
                                    int method )

{
  // follow same argument scheme as InverseCompton.cc for compatibility

  if(debug==1) cout<<">> Bremsstrahlung:bremss_emiss()"<< " method="<<method<<endl;

  // assuming electron energy scale is logarithmic

 
  emiss.resize(Egamma.size());
  emiss=0.0;

 if(method==1)
 {

  double bremss_sigma;

  

  for (int iEgamma=0;iEgamma<Egamma.size();iEgamma++)
  for (int iEe    =0;iEe    <Ee    .size();iEe    ++)
  {
    bremss_sigma      = dSdE (Ee[iEe], Egamma[iEgamma],IZ1,Ne);

    emiss[iEgamma]   += bremss_sigma * fluxe[iEe] 
                                     *    Ee[iEe] ;

    if(debug==1) cout<<" bremss_emiss: method="<<method<<" Egamma="<< Egamma[iEgamma]<<" Ee="<<Ee[iEe]
                     <<" bremss_sigma= "<<bremss_sigma<<" emiss="<<emiss[iEgamma]<<endl;
  }


  emiss *= log(Ee[1]/Ee[0]);// NB assumes at least 2 values in spectra  !
  emiss *= 1.0e-24; // cross-section in barn MeV-1

 }// method==1

 return;
}

///////////////////////////////////////////////////////////////////////////

 valarray<double >Bremsstrahlung:: bremss_emiss(valarray<double> Ee,   valarray<double> Egamma,
                                    valarray<double> fluxe,
                                    int IZ1, int Ne,
                                    int method )

{
  // same as previous function but returns the array instead of using argument

  if(debug==1) cout<<">> Bremsstrahlung:valarray<double> bremss_emiss()"<< " method="<<method<<endl;

  // assuming electron energy scale is logarithmic

 
  valarray<double> emiss;
  emiss.resize(Egamma.size());
  emiss=0.0;

 if(method==1)
 {

  double bremss_sigma;

  
  for (int iEgamma=0;iEgamma<Egamma.size();iEgamma++)
  for (int iEe    =0;iEe    <Ee    .size();iEe    ++)
  {
    bremss_sigma      = dSdE (Ee[iEe], Egamma[iEgamma],IZ1,Ne);

    emiss[iEgamma]   += bremss_sigma * fluxe[iEe] 
                                     *    Ee[iEe] ;

    if(debug==1) cout<<" bremss_emiss valarray<double> : method="<<method<<" Egamma="<< Egamma[iEgamma]
                     <<" Ee="<<Ee[iEe] <<" fluxe="<<fluxe[iEe]
                     <<" bremss_sigma= "<<bremss_sigma<<" emiss="<<emiss[iEgamma]<<endl;
  }


  emiss *= log(Ee[1]/Ee[0]);// NB assumes at least 2 values in spectra  !
  emiss *= 1.0e-24; // cross-section in barn MeV-1

 }// method==1

   if(debug==1)
  for (int iEgamma=0;iEgamma<Egamma.size();iEgamma++)
   cout<<" bremss_emiss valarray<double> : method="<<method<<" Egamma="<< Egamma[iEgamma]
       <<" emiss="<<emiss[iEgamma]<<endl;


 return emiss;
}

///////////////////////////////////////////////////////////////////

double  Bremsstrahlung::emissivity_band
(valarray<double> Ee, valarray<double> fluxe,
                         int IZ1, int Ne,
                         double E_gamma_band_min, double E_gamma_band_max,
                         double E_gamma_factor,
                         int method )
{

  if(debug==1)cout<<"Bremsstrahlung::emissivity_band"<<endl;

  valarray<double> E_gamma_work;
  

  int nE_gamma_work=log(E_gamma_band_max/E_gamma_band_min)/log(E_gamma_factor)+1.001; 

  E_gamma_work.resize(nE_gamma_work);

  for(int iE_gamma_work =0;iE_gamma_work <nE_gamma_work  ;iE_gamma_work++)
        E_gamma_work[iE_gamma_work] = E_gamma_band_min * pow(E_gamma_factor,iE_gamma_work);

  valarray<double> bremss_emiss_;
  bremss_emiss_.resize(nE_gamma_work);

  bremss_emiss_ =   bremss_emiss( Ee,    E_gamma_work,       fluxe, IZ1,  Ne,    method );

 double emiss_band=0.;

 for(int iE_gamma_work =0;iE_gamma_work <nE_gamma_work  ;iE_gamma_work++ )
    {
      emiss_band+= bremss_emiss_  [iE_gamma_work ]*E_gamma_work[iE_gamma_work ];
      
      if(debug==1)
      cout<<"emissivity_band single, cumulative sum: "
          <<" E_gamma_band_min="<< E_gamma_band_min
          <<" E_gamma_band_max="<< E_gamma_band_max
          <<" E_gamma= "<<E_gamma_work[iE_gamma_work]<<" MeV "
	  <<" bremss_emiss_= "<<  bremss_emiss_ [iE_gamma_work]
                                                   <<endl;
    }

  emiss_band*=log(E_gamma_work[1]/E_gamma_work[0]);

  if(debug==1)
  cout<<"emissivity_band single value:"                
             <<" E_gamma_band_min="<< E_gamma_band_min
             <<" E_gamma_band_max="<< E_gamma_band_max
             <<" emiss_band="<<emiss_band<<endl;

  return emiss_band;
};

///////////////////////////////////////////////////////////////////

valarray<double>  Bremsstrahlung::emissivity_band
                        (valarray<double> Ee, valarray<double> fluxe,
                         int IZ1, int Ne,
                         valarray<double> E_gamma_band_min,
                         valarray<double> E_gamma_band_max, 
                         double E_gamma_factor,
                         int method )
{

  if(debug==1)cout<<"Bremsstrahlung::valarray<double> emissivity_band "<<endl;

  int nE_gamma_band=E_gamma_band_min.size();
  valarray<double> emiss_band;
  emiss_band.resize(nE_gamma_band);
  emiss_band=0.;
  
  for (int iE_gamma_band=0; iE_gamma_band< nE_gamma_band; iE_gamma_band++)
  {
    
   emiss_band[iE_gamma_band] = emissivity_band
                               ( Ee, fluxe,
                                 IZ1,  Ne,
                                 E_gamma_band_min[iE_gamma_band], E_gamma_band_max[iE_gamma_band],E_gamma_factor,
                                 method);
    
  if(debug==1)
         cout<<"emissivity_band valarray<double>:"
             <<" E_gamma_band_min="<< E_gamma_band_min[iE_gamma_band]
             <<" E_gamma_band_max="<< E_gamma_band_max[iE_gamma_band]
             <<" emiss_band="<<emiss_band[iE_gamma_band]<<endl;

  }//  iE_gamma_band

  return emiss_band;
};

/////////////////////////////////////////////////////////////

valarray<double>  Bremsstrahlung:: emissivity_band 
                                  ( valarray<double> fluxe)
{

  int debug=1;
  if(debug==1)cout<<"Bremsstrahlung::emissivity_band using matrices"<<endl;

  int nE_gamma_band=gamma_band_matrix[0].size();
  valarray<double> emiss_band;
  emiss_band.resize(nE_gamma_band);
  emiss_band=0.;

  for (int iE_gamma_band=0; iE_gamma_band< nE_gamma_band; iE_gamma_band++)
  {
   for(int iEe=0;iEe<fluxe.size();  iEe++)
   {
    emiss_band[iE_gamma_band] += 
      gamma_band_matrix[iEe][ iE_gamma_band] * fluxe[iEe] * Ee_grid[iEe];

   }//iEe
  }// iE_gamma_band

  emiss_band *=log(Ee_grid[1]/Ee_grid[0]); // log integration
  return emiss_band;
}

///////////////////////////////////////////////////////////

int  Bremsstrahlung::gen_gamma_band_matrices(valarray<double> Ee, 
                                             valarray<int> IZ1_list, valarray<int> Ne_list, 
                                             valarray<double> target_abundances,
                                             valarray<double> E_gamma_band_min,
                                             valarray<double> E_gamma_band_max, 
                                             double E_gamma_factor,
                                             int method )

{
  int debug=1;
  if(debug==1)cout<<"Bremsstrahlung::gen_gamma_band_matrices"<<endl;

  gamma_band_matrix .resize(Ee.size());
  int nE_gamma_band=E_gamma_band_min.size();

  for(int iEe=0;iEe<Ee.size();iEe++)
  {
   gamma_band_matrix [iEe].resize(nE_gamma_band);
  }

  
  if(debug>=1)
  for(int iE_gamma_band=0; iE_gamma_band< nE_gamma_band; iE_gamma_band++)
         cout<<"gen_gamma_band_matrices: input"
             <<" E_gamma_band_min="<< E_gamma_band_min[iE_gamma_band]
             <<" E_gamma_band_max="<< E_gamma_band_max[iE_gamma_band]
             <<endl;


    
 for(int iEe=0;iEe<Ee.size();iEe++)
 {
  for(int iE_gamma_band=0; iE_gamma_band< nE_gamma_band; iE_gamma_band++)
  {
   int nE_gamma_work=log(E_gamma_band_max[iE_gamma_band]/E_gamma_band_min[iE_gamma_band])/log(E_gamma_factor)+1.001; 
   valarray<double> E_gamma_work;   E_gamma_work.resize(nE_gamma_work);

   for(int iE_gamma_work =0;iE_gamma_work <nE_gamma_work  ;iE_gamma_work++)
        E_gamma_work[iE_gamma_work] = E_gamma_band_min[iE_gamma_band] * pow(E_gamma_factor,iE_gamma_work);

   double emiss_band=0.;

   for(int iE_gamma_work =0;iE_gamma_work <nE_gamma_work  ;iE_gamma_work++ )
   {
     double bremss_sigma=0;
     for(int i=0;i<IZ1_list.size();i++)
            bremss_sigma += dSdE (Ee[iEe], E_gamma_work[iE_gamma_work],IZ1_list[i],Ne_list[i]) * target_abundances[i];

     emiss_band+=  bremss_sigma * E_gamma_work[iE_gamma_work ];
      
      if(debug>=3)
      cout<<"gen_gamma_band_matrices single, cumulative sum: "
          <<" Ee="<<Ee[iEe]
          <<" E_gamma_band_min="<< E_gamma_band_min[iE_gamma_band]
          <<" E_gamma_band_max="<< E_gamma_band_max[iE_gamma_band]
          <<" E_gamma= "<<E_gamma_work[iE_gamma_work]<<" MeV "
	  <<" emiss_band= "<<  emiss_band
          <<endl;
    }

   emiss_band *= log(E_gamma_work[1]/E_gamma_work[0]);  // log integration 
   emiss_band *= 1.0e-24; // cross-section in barn MeV-1, convert to cm^2 MeV-1

  if(debug>=2)
  cout<<"gen_gamma_band_matrices single value:"                
             <<" E_gamma_band_min="<< E_gamma_band_min[iE_gamma_band]
             <<" E_gamma_band_max="<< E_gamma_band_max[iE_gamma_band]
             <<" emiss_band="<<emiss_band<<endl;

  gamma_band_matrix [iEe][iE_gamma_band] = emiss_band;

 
   }//iE_gamma_band
  }//iEe

 if(debug>=1)
             for(int iEe=0;iEe<Ee.size();iEe++)
             for(int iE_gamma_band=0; iE_gamma_band< nE_gamma_band; iE_gamma_band++)
             cout<<"gen_gamma_band_matrices result:" 
             <<" E_gamma_band_min="<< E_gamma_band_min[iE_gamma_band]
             <<" E_gamma_band_max="<< E_gamma_band_max[iE_gamma_band]
             <<" gamma_band_matrix= "<<  gamma_band_matrix [iEe][iE_gamma_band]
             <<endl;

 // store the electron energy grid in the class for use in integration
 Ee_grid.resize(Ee.size());
 Ee_grid = Ee;


  return 0;
}






////////////////////////////////////////////////// test program ////////////////////

void Bremsstrahlung::test(int test_dSdE, int test_emiss)
{
  cout<< " Bremsstrahlung::test"<<endl;

  string name_="test name";
  valarray<double> parameters_;
  int debug_=1;

  init(name_,parameters_,debug_);

  double Ee,Egam; int IZ1, Ne1; double dSdK;

 

  if( test_dSdE==1)
  {

  Egam=   10.; // gamma-ray energy, MeV
  Ee  = 1000.; // electron energy, MeV 
  IZ1 = 1;     // nucleus charge
  Ne1 = 0;     // 0=unshielded, 1,2=shielded 1,2 electron atoms

  /*
  Egam=350.343; 
  Ee  =352.224;
  double dSdE_=dSdE( Ee, Egam, IZ1,  Ne1 );
  cout<<  " dSdE debug test case Ee="<<Ee<<" Egam="<<Egam<< " IZ1="<<IZ1<<"  Ne1="<<Ne1<<" dSdE="<<dSdE_<<endl;

  Ee  = 369.835;
  dSdE_=dSdE( Ee, Egam, IZ1,  Ne1 );
  cout<<  " dSdE debug test case Ee="<<Ee<<" Egam="<<Egam<< " IZ1="<<IZ1<<"  Ne1="<<Ne1<<" dSdE="<<dSdE_<<endl;
  */


 // test each part of dSdE
  Ee=0.4;     Egam=0.01;  IZ1=1; Ne1=0; double dSdE_=dSdE( Ee, Egam, IZ1,  Ne1 );
  cout<<  "dSdE result part 0: Ee="<<Ee<<" Egam="<<Egam<< " IZ1="<<IZ1<<"  Ne1="<<Ne1<<" dSdE="<<dSdE_<<endl;
  Ee=0.523284;Egam=0.01;  IZ1=1; Ne1=0;        dSdE_=dSdE( Ee, Egam, IZ1,  Ne1 );
  cout<<  "dSdE result part 1: Ee="<<Ee<<" Egam="<<Egam<< " IZ1="<<IZ1<<"  Ne1="<<Ne1<<" dSdE="<<dSdE_<<endl;
  Ee=0.583811; Egam=0.01; IZ1=2; Ne1=0;        dSdE_=dSdE( Ee, Egam, IZ1,  Ne1 );
  cout<<  "dSdE result part 2: Ee="<<Ee<<" Egam="<<Egam<< " IZ1="<<IZ1<<"  Ne1="<<Ne1<<" dSdE="<<dSdE_<<endl;
  Ee=2.52066;  Egam=0.01; IZ1=3; Ne1=3;        dSdE_=dSdE( Ee, Egam, IZ1,  Ne1 );
  cout<<  "dSdE result part 3: Ee="<<Ee<<" Egam="<<Egam<< " IZ1="<<IZ1<<"  Ne1="<<Ne1<<" dSdE="<<dSdE_<<endl;
  Ee=2.52066;  Egam=0.01; IZ1=2;Ne1=2;        dSdE_=dSdE( Ee, Egam, IZ1,  Ne1 );
  cout<<  "dSdE result part 4: Ee="<<Ee<<" Egam="<<Egam<< " IZ1="<<IZ1<<"  Ne1="<<Ne1<<" dSdE="<<dSdE_<<endl;
  Ee=2.52066  ;Egam=0.01; IZ1=2;Ne1=0;        dSdE_=dSdE( Ee, Egam, IZ1,  Ne1 );
  cout<<  "dSdE result part 5: Ee="<<Ee<<" Egam="<<Egam<< " IZ1="<<IZ1<<"  Ne1="<<Ne1<<" dSdE="<<dSdE_<<endl;
  Ee=2.52066;  Egam=0.01; IZ1=2; Ne1=1;       dSdE_=dSdE( Ee, Egam, IZ1,  Ne1 );
  cout<<  "dSdE result part 6: Ee="<<Ee<<" Egam="<<Egam<< " IZ1="<<IZ1<<"  Ne1="<<Ne1<<" dSdE="<<dSdE_<<endl;
  Ee=2.52066;  Egam=0.01; IZ1=1; Ne1=1;       dSdE_=dSdE( Ee, Egam, IZ1,  Ne1 );
  cout<<  "dSdE result part 7: Ee="<<Ee<<" Egam="<<Egam<< " IZ1="<<IZ1<<"  Ne1="<<Ne1<<" dSdE="<<dSdE_<<endl;


  for(Ee  =4.e-1;Ee  <1e6;Ee  *=1.01)// NB total electron energy
  for(Egam=1.e-2;Egam<1e6;Egam*=2.  )
  {

  for(IZ1=1;IZ1<=4;IZ1++)
  {
   for (Ne1=0;Ne1<=4;Ne1++)
   {
      double dSdE_=dSdE( Ee, Egam, IZ1,  Ne1 );

      cout<<  "dSdE result: Ee="<<Ee<<" Egam="<<Egam<< " IZ1="<<IZ1<<"  Ne1="<<Ne1<<" dSdE="<<dSdE_<<endl;
     }
   }
  }

 

}// test_dSdE

  ///////////////////////////////////////////////////////

 

 if(test_emiss==1)
  {

  cout<<" testing bremsstrahlung emissivity"<<endl;



  double  Eemin       =10.0; // MeV
  double  Eemax       = 1.0e6; 
  double  Eefactor    = 1.05;

  double  Egammamin   = 1.0; // MeV
  double  Egammamax   = 1.0e6; 
  double  Egammafactor=pow(10.,0.1);

  int nEe    =log(    Eemax/    Eemin)/log(Eefactor)    +1.0001;// NB must be >1 for correct run 
  int nEgamma=log(Egammamax/Egammamin)/log(Egammafactor)+1.0001; 


  double ge=3.0; // electron spectral index
 
  valarray<double> Ee_grid      (nEe);
  valarray<double> fluxe        (nEe);
  valarray<double> Egamma_grid  (nEgamma);
  valarray<double> emiss10      (nEgamma);
  valarray<double> emiss11      (nEgamma);
  valarray<double> emiss20      (nEgamma);
  valarray<double> emiss21      (nEgamma);
  valarray<double> emiss22      (nEgamma);

  int IZ1,Ne;
  IZ1=1; // target nucleus charge
  Ne=0;  // 0=unscreened, 1=1-electron atom, 2=2-electron atom

  int method = 1;
  
  debug_=0;
  init(name_,parameters_,debug_);
  

  for (int iEgamma=0;iEgamma<Egamma_grid.size();iEgamma++) Egamma_grid [iEgamma]= Egammamin*pow(Egammafactor,iEgamma);
  for (int iEe    =0;iEe    <Ee_grid    .size();iEe++)         Ee_grid [iEe]    = Eemin    *pow(Eefactor,    iEe    );
   

  // Ackermann etal Fermi-LAT arXiv:1109.052 electrons fit, converting from GeV m^2 to MeV cm^2 

  ge=3.19;
  
  fluxe     = 2.07e-9* pow(Ee_grid/2.e4, -ge   );

  //  ge=2.0; //test, effective break for low energy electrons
  //  fluxe     = 2.07e-9* pow(Ee_grid/2.e4, -ge   );

  for (int iEe    =0;iEe    <Ee_grid    .size();iEe++) cout<<"electron spectrum: Ee= "<< Ee_grid    [iEe]<<" fluxe="<<fluxe[iEe]<<endl;


  // differential emissivity spectrum

  for(method=1;method<=1;method++)
  {

  IZ1=1;  // target nucleus charge
  Ne =0;  // 0=unscreened, 1=1-electron atom, 2=2-electron atom

   bremss_emiss(Ee_grid,  Egamma_grid, fluxe, IZ1,Ne, emiss10, method );

   IZ1=1;
   Ne =1;
   bremss_emiss(Ee_grid,  Egamma_grid, fluxe, IZ1,Ne, emiss11, method );

   IZ1=2;// should be  >=Ne 
   Ne =0;
   bremss_emiss(Ee_grid,  Egamma_grid, fluxe, IZ1,Ne, emiss20, method );

   IZ1=2;// should be  >=Ne 
   Ne =1;
   bremss_emiss(Ee_grid,  Egamma_grid, fluxe, IZ1,Ne, emiss21, method );

   IZ1=2;// should be  >=Ne 
   Ne =2;
   bremss_emiss(Ee_grid,  Egamma_grid, fluxe, IZ1,Ne, emiss22, method );

 for (int iEgamma=0;iEgamma<Egamma_grid.size();iEgamma++)
       cout<<" bremss_emiss argument version, method="<<method<<" Egamma= "<< Egamma_grid[iEgamma]
           <<" bremss emiss10="<<emiss10[iEgamma] 
           <<" bremss emiss11="<<emiss11[iEgamma] 
           <<" bremss emiss20="<<emiss20[iEgamma]
           <<" bremss emiss21="<<emiss21[iEgamma]
           <<" bremss emiss22="<<emiss22[iEgamma]
           <<endl;


  }//method
  

  for(method=1;method<=1;method++)
  {

  IZ1=1;  // target nucleus charge
  Ne =0;  // 0=unscreened, 1=1-electron atom, 2=2-electron atom

   emiss10=bremss_emiss(Ee_grid,  Egamma_grid, fluxe, IZ1,Ne, method );

   IZ1=1;
   Ne =1;
   emiss11=bremss_emiss(Ee_grid,  Egamma_grid, fluxe, IZ1,Ne, method );

   IZ1=2;// should be  >=Ne 
   Ne =0;
   emiss20=bremss_emiss(Ee_grid,  Egamma_grid, fluxe, IZ1,Ne,  method );

   IZ1=2;// should be  >=Ne 
   Ne =1;
   emiss21=bremss_emiss(Ee_grid,  Egamma_grid, fluxe, IZ1,Ne, method );

   IZ1=2;// should be  >=Ne 
   Ne =2;
   emiss22=bremss_emiss(Ee_grid,  Egamma_grid, fluxe, IZ1,Ne, method );

 for (int iEgamma=0;iEgamma<Egamma_grid.size();iEgamma++)
       cout<<" bremss_emiss valarray version, method="<<method<<" Egamma= "<< Egamma_grid[iEgamma]
           <<" bremss emiss10="<<emiss10[iEgamma] 
           <<" bremss emiss11="<<emiss11[iEgamma] 
           <<" bremss emiss20="<<emiss20[iEgamma]
           <<" bremss emiss21="<<emiss21[iEgamma]
           <<" bremss emiss22="<<emiss22[iEgamma]
           <<endl;


  }//method




  // emissivity in energy bands

  double E_gamma_band_min   =  100.;// MeV
  double E_gamma_band_max   =10000.;
  double E_gamma_band_factor= 1.02;// for integration over gamma-ray band

  IZ1=1;  // target nucleus charge
  Ne =0;  // 0=unscreened, 1=1-electron atom, 2=2-electron atom
  Ne =1;

  /*
  IZ1=2;  // target nucleus charge
  Ne =1;  // 0=unscreened, 1=1-electron atom, 2=2-electron atom
  */

  method=1;

  double emiss_band=emissivity_band( Ee_grid, fluxe,
                                     IZ1,  Ne,
                                     E_gamma_band_min, E_gamma_band_max,
                                     E_gamma_band_factor,
		                     method );

      cout<<"Ackermann etal electrons "                 
          <<" E_gamma="<< E_gamma_band_min<<" -"<< E_gamma_band_max
	  <<" MeV predicted emiss="  <<emiss_band   <<endl;


  int             nE_gamma_band;
  valarray<double> E_gamma_band_min_;
  valarray<double> E_gamma_band_max_;
  valarray<double> emiss_band_HI; // atomic H
  valarray<double> emiss_band_He; // atomic He
  valarray<double> emiss_band_;


 // for cf values in Fermi 2nd quadrant paper http://adsabs.harvard.edu/abs/2010ApJ...710..133A  Table 1 qHI,1=local
   // qHI=0.2-0.4 GeV: 0.584,  0.4-0.6 GeV:0.224, 0.6-1.0 GeV: 0.168,  1-2 GeV : 0.110, 2-10 GeV :0.048e-26 sr sr-1

  nE_gamma_band=35;
  E_gamma_band_min_.resize(nE_gamma_band);
  E_gamma_band_max_.resize(nE_gamma_band);

  
  E_gamma_band_min_[ 0]= 0.050; E_gamma_band_max_[ 0]= 0.100;//GeV
  E_gamma_band_min_[ 1]= 0.100; E_gamma_band_max_[ 1]= 0.200;

  // 2nd quadrant paper bands
  E_gamma_band_min_[ 2]= 0.200; E_gamma_band_max_[ 2]= 0.400;
  E_gamma_band_min_[ 3]= 0.400; E_gamma_band_max_[ 3]= 0.600;
  E_gamma_band_min_[ 4]= 0.600; E_gamma_band_max_[ 4]= 1.000;
  E_gamma_band_min_[ 5]= 1.000; E_gamma_band_max_[ 5]= 2.000;
  E_gamma_band_min_[ 6]= 2.000; E_gamma_band_max_[ 6]=10.000;

  // another band
  E_gamma_band_min_[ 7]=10.000; E_gamma_band_max_[ 7]=20.000;

  // Abdo etal 2009 HI emissivity paper bands
  E_gamma_band_min_[ 8]= 0.100; E_gamma_band_max_[ 8]= 0.140;//GeV
  E_gamma_band_min_[ 9]= 0.140; E_gamma_band_max_[ 9]= 0.200;
  E_gamma_band_min_[10]= 0.200; E_gamma_band_max_[10]= 0.280;
  E_gamma_band_min_[11]= 0.280; E_gamma_band_max_[11]= 0.400;
  E_gamma_band_min_[12]= 0.400; E_gamma_band_max_[12]= 0.560;
  E_gamma_band_min_[13]= 0.560; E_gamma_band_max_[13]= 0.800;
  E_gamma_band_min_[14]= 0.800; E_gamma_band_max_[14]= 1.130;
  E_gamma_band_min_[15]= 1.130; E_gamma_band_max_[15]= 1.600;
  E_gamma_band_min_[16]= 1.600; E_gamma_band_max_[16]= 2.260;
  E_gamma_band_min_[17]= 2.260; E_gamma_band_max_[17]= 3.200;
  E_gamma_band_min_[18]= 3.200; E_gamma_band_max_[18]= 4.530;
  E_gamma_band_min_[19]= 4.530; E_gamma_band_max_[19]= 6.400;
  E_gamma_band_min_[20]= 6.400; E_gamma_band_max_[20]= 9.050;

  // Ackermann etal 2011 3rd quadrant paper
  E_gamma_band_min_[21]= 0.100; E_gamma_band_max_[21]= 0.140;//GeV
  E_gamma_band_min_[22]= 0.140; E_gamma_band_max_[22]= 0.200;
  E_gamma_band_min_[23]= 0.200; E_gamma_band_max_[23]= 0.280;
  E_gamma_band_min_[24]= 0.280; E_gamma_band_max_[24]= 0.400;
  E_gamma_band_min_[25]= 0.400; E_gamma_band_max_[25]= 0.560;
  E_gamma_band_min_[26]= 0.560; E_gamma_band_max_[26]= 0.800;
  E_gamma_band_min_[27]= 0.800; E_gamma_band_max_[27]= 1.130;
  E_gamma_band_min_[28]= 1.130; E_gamma_band_max_[28]= 1.600;
  E_gamma_band_min_[29]= 1.600; E_gamma_band_max_[29]= 2.260;
  E_gamma_band_min_[30]= 2.260; E_gamma_band_max_[30]= 3.200;
  E_gamma_band_min_[31]= 3.200; E_gamma_band_max_[31]= 4.530;
  E_gamma_band_min_[32]= 4.530; E_gamma_band_max_[32]= 6.400;
  E_gamma_band_min_[33]= 6.400; E_gamma_band_max_[33]= 9.050;
  E_gamma_band_min_[34]= 9.050; E_gamma_band_max_[34]=25.600;

  E_gamma_band_min_    *=1000.;// the values above are in GeV, Bremsstrahlung works in MeV
  E_gamma_band_max_    *=1000.;// the values above are in GeV, Bremsstrahlung works in MeV


  valarray<double> emiss_measured,emiss_measured_from_table;
  emiss_measured           .resize(nE_gamma_band);
  emiss_measured = 0.0;
  emiss_measured[ 2]=0.584e-26;
  emiss_measured[ 3]=0.224e-26;
  emiss_measured[ 4]=0.168e-26;
  emiss_measured[ 5]=0.110e-26;
  emiss_measured[ 6]=0.048e-26;
  for(int i=2;i<= 6;i++) 
   cout<<"Abdo etal 2010 ApJ 710, 133 2nd quadrant Goulds Belt measured emissivity: E_gamma="<<E_gamma_band_min_[i]
        <<" - "<<E_gamma_band_max_[i]
        <<" GeV =" << emiss_measured[i]<<endl;
  


  // Table 1 of Abdo etal 703, 1249, 2009 (HI emissivity paper)  MeV^2 s-1 sr-1 MeV-1
  emiss_measured[ 8] = 1.04e-24;
  emiss_measured[ 9] = 1.67e-24;
  emiss_measured[10] = 1.91e-24;
  emiss_measured[11] = 2.11e-24;
  emiss_measured[12] = 2.33e-24;
  emiss_measured[13] = 2.20e-24;
  emiss_measured[14] = 2.17e-24;
  emiss_measured[15] = 1.88e-24;
  emiss_measured[16] = 1.72e-24;
  emiss_measured[17] = 1.15e-24;
  emiss_measured[18] = 1.10e-24;
  emiss_measured[19] = 1.12e-24;
  emiss_measured[20] = 0.71e-24;
  

  // convert from MeV^2 s-1 sr-1 MeV-1 to  s-1 sr-1 using geometric mean for E^2
  for(int i=8;i<=20;i++)
  {
    //emiss_measured[i]*=1.0e-3; // MeV required so don't convert to GeV here
   emiss_measured[i]*=(E_gamma_band_max_[i]- E_gamma_band_min_[i]);
   emiss_measured[i]/=(E_gamma_band_max_[i]* E_gamma_band_min_[i]);
   cout<<"Abdo etal 2009 703, 1249, 2009 (HI emissivity paper) measured emissivity: E_gamma="
        <<E_gamma_band_min_[i]<<" - "<<E_gamma_band_max_[i]
        <<" MeV =" << emiss_measured[i]<<endl;
  }

  // Table 2 of Ackermann etal ApJ 726, 81, 2011   (3rd quad paper)  MeV^2 s-1 sr-1 MeV-1
  // Ts=250K, qHI,1 = local
  emiss_measured[21] = 1.35e-24;
  emiss_measured[22] = 1.56e-24;
  emiss_measured[23] = 1.82e-24;
  emiss_measured[24] = 2.00e-24;
  emiss_measured[25] = 1.95e-24;
  emiss_measured[26] = 2.10e-24;
  emiss_measured[27] = 1.90e-24;
  emiss_measured[28] = 1.79e-24;
  emiss_measured[29] = 1.74e-24;
  emiss_measured[30] = 1.37e-24;
  emiss_measured[31] = 0.84e-24;
  emiss_measured[32] = 0.65e-24;
  emiss_measured[33] = 0.54e-24;
  emiss_measured[34] = 0.38e-24;

  // convert from MeV^2 s-1 sr-1 MeV-1 to  s-1 sr-1 using geometric mean for E^2
  for(int i=21;i<=34;i++)
  {
    // emiss_measured[i]*=1.0e-3; //  MeV required so don't convert to GeV here
   emiss_measured[i]*=(E_gamma_band_max_[i]- E_gamma_band_min_[i]);
   emiss_measured[i]/=(E_gamma_band_max_[i]* E_gamma_band_min_[i]);
   cout<<"Ackermann etal ApJ 726, 81, 2011   (3rd quad paper) measured emissivity: E_gamma="
        <<E_gamma_band_min_[i]<<" - "<<E_gamma_band_max_[i]
        <<" MeV =" << emiss_measured[i]<<endl;
  }


  emiss_band_HI.resize(nE_gamma_band);
  emiss_band_He.resize(nE_gamma_band);
  emiss_band_  .resize(nE_gamma_band);


  emiss_band_=emissivity_band( Ee_grid, fluxe,
                                     IZ1,  Ne,
                                     E_gamma_band_min_, E_gamma_band_max_,
                                     E_gamma_band_factor,
		                     method );

 for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
             cout<<" Bremsstrahlung, Ackermann etal  electrons, output from multiple bands version:  "                 
             <<" E_gamma="<< E_gamma_band_min_[iE_gamma_band]<<" -"<< E_gamma_band_max_[iE_gamma_band]
             <<" MeV predicted emiss="  <<emiss_band_[iE_gamma_band] 
	     <<" measured = "           <<emiss_measured[iE_gamma_band]
	     <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;


  CR_Spectrum cr_spectrum;
  string cr_spectrum_name="Fermi2011_electrons_positrons";

  valarray<double>cr_spectrum_parameters;
  cr_spectrum_parameters.resize(10);

  double Phi0=0;
  double  ge1=2.;//electron index below break


  double He_fraction=0.1; // He/H by number

  for(ge1=1.0;ge1<=3.2;ge1+=0.2)
  {
  cr_spectrum_parameters[0]=Phi0;
  cr_spectrum_parameters[1]=3.e3;// Ee_break
  cr_spectrum_parameters[2]=ge1; // ge1=electron index below break
  cr_spectrum_parameters[3]=5.e3;// Ep_break
  cr_spectrum_parameters[4]=2.0; // gp1=positron index below break

  cr_spectrum_parameters[5]=1.0;//  electrons factor
  cr_spectrum_parameters[6]=1.0;//  positrons factor

  int cr_spectrum_debug=1;
  cr_spectrum.init(cr_spectrum_name,cr_spectrum_parameters,cr_spectrum_debug);
  

  double  r=1.0;// not used for "Fermi2011_electrons_positrons" but argument required
  fluxe=cr_spectrum.fluxe(Ee_grid,r);

  for(int i=0;i<Ee_grid.size();i++) cout<<  " cr_spectrum: Ee= "<<  Ee_grid[i]<<" fluxe="<<fluxe[i]<<endl;

  E_gamma_band_factor= 1.02;// for integration over gamma-ray band

  IZ1=1;Ne=1;
  emiss_band_HI=emissivity_band( Ee_grid, fluxe,
                                     IZ1,  Ne,
                                     E_gamma_band_min_, E_gamma_band_max_,
                                     E_gamma_band_factor,
		                     method );

  IZ1=2;Ne=2;
  emiss_band_He=emissivity_band( Ee_grid, fluxe,
                                     IZ1,  Ne,
                                     E_gamma_band_min_, E_gamma_band_max_,
                                     E_gamma_band_factor,
		                     method );


  emiss_band_He *= He_fraction;
  emiss_band_ = emiss_band_HI + emiss_band_He;


  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
         cout<<" Bremsstrahlung: "
             <<" He fraction="<<He_fraction
	     <<" CR spectrum="<<cr_spectrum_name
             <<" electron index below break ="<<ge1 
             <<" E_gamma="<< E_gamma_band_min_[iE_gamma_band]<<" -"<< E_gamma_band_max_[iE_gamma_band]
             <<" MeV predicted emiss HI="  <<emiss_band_HI[iE_gamma_band] 
             <<                    " He="  <<emiss_band_He[iE_gamma_band] 
             <<                 " total="  <<emiss_band_  [iE_gamma_band] 
	     <<" measured = "             <<emiss_measured[iE_gamma_band]
	     <<" predicted/measured = "   <<emiss_band_   [iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;


  // test generation of matrices
 
 int n_target=1;
 valarray<int> IZ1_list;
 valarray<int>  Ne_list;
 valarray<double> target_abundances;

 IZ1_list.resize(n_target);
  Ne_list.resize(n_target);
 target_abundances.resize(n_target);


 
 method=1;

 valarray<double>emiss_band_from_matrices_HI;
 valarray<double>emiss_band_from_matrices_He;
 valarray<double>emiss_band_from_matrices;
 valarray<double>emiss_band_from_matrices_4;
 emiss_band_from_matrices_HI.resize(nE_gamma_band);
 emiss_band_from_matrices_He.resize(nE_gamma_band);
 emiss_band_from_matrices   .resize(nE_gamma_band);
 emiss_band_from_matrices_4 .resize(nE_gamma_band);

 IZ1_list[0]=1; Ne_list[0]=1;   target_abundances[0]=1.0;

 gen_gamma_band_matrices( Ee_grid, 
                          IZ1_list,  Ne_list,  target_abundances,
			  E_gamma_band_min_,E_gamma_band_max_,E_gamma_band_factor,                                                                    method );

 emiss_band_from_matrices_HI = emissivity_band(fluxe);


 IZ1_list[0]=2; Ne_list[0]=2;   target_abundances[0]=1.0;

 gen_gamma_band_matrices( Ee_grid, 
                          IZ1_list,  Ne_list,  target_abundances,
			  E_gamma_band_min_,E_gamma_band_max_,E_gamma_band_factor,                                                                    method );

 emiss_band_from_matrices_He = emissivity_band(fluxe);

 n_target=2;
 IZ1_list.resize(n_target);
  Ne_list.resize(n_target);
 target_abundances.resize(n_target);
 IZ1_list[0]=1; Ne_list[0]=1; target_abundances[0]=1.0;
 IZ1_list[1]=2; Ne_list[1]=2; target_abundances[1]=0.1;// He fraction by number



 gen_gamma_band_matrices( Ee_grid, 
                          IZ1_list,  Ne_list,  target_abundances,
			  E_gamma_band_min_,E_gamma_band_max_,E_gamma_band_factor,                                                                    method );

 emiss_band_from_matrices    = emissivity_band(fluxe);


 n_target=4;
 IZ1_list.resize(n_target);
  Ne_list.resize(n_target);
 target_abundances.resize(n_target);
 IZ1_list[0]=1; Ne_list[0]=1; target_abundances[0]=1.0; // HI
 IZ1_list[1]=2; Ne_list[1]=2; target_abundances[1]=0.1 ;// He fraction by number
 IZ1_list[2]=1; Ne_list[2]=0; target_abundances[2]=0.1 ;// ionized H just test value
 IZ1_list[3]=2; Ne_list[3]=1; target_abundances[3]=0.05;// ionized He just test value

 gen_gamma_band_matrices( Ee_grid, 
                          IZ1_list,  Ne_list,  target_abundances,
			  E_gamma_band_min_,E_gamma_band_max_,E_gamma_band_factor,                                                                    method );

 emiss_band_from_matrices_4  = emissivity_band(fluxe);

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
         cout<<" emiss using matrices: "
             <<" He fraction="<<He_fraction
	     <<" CR spectrum="<<cr_spectrum_name
             <<" electron index below break ="<<ge1 
             <<" E_gamma="<< E_gamma_band_min_[iE_gamma_band]<<" -"<< E_gamma_band_max_[iE_gamma_band]

	     <<" from matrices HI="    <<emiss_band_from_matrices_HI[iE_gamma_band] 
	     <<              " He="    <<emiss_band_from_matrices_He[iE_gamma_band] 
	     <<              " HI+He=" <<emiss_band_from_matrices   [iE_gamma_band]
	     <<              " HI+He+HII+HeII =" <<emiss_band_from_matrices_4 [iE_gamma_band]
	     <<" direct HI="  <<emiss_band_HI[iE_gamma_band] 
	     <<" He="         <<emiss_band_He[iE_gamma_band] 
             <<" HI+He="      <<emiss_band_  [iE_gamma_band] 
             <<" ratio of HI+He ="<<emiss_band_from_matrices[iE_gamma_band]/emiss_band_[iE_gamma_band] 
	     <<" measured = "               <<emiss_measured[iE_gamma_band]
	     <<" predicted/measured = "     <<emiss_band_   [iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;



  }// ge1



  }// if test_emiss==1



  return ;
};


///////////////////////////////////////////////////////////////////////////////////////////

/* Original fortran documentation

      subroutine bremss_spec(Egam,E0,IZ1,Ne1,dSdK) 
c***********************************************************************
c          *** I.Moskalenko (MPE,Garching) *** version of 14.05.1998 ***
c PURPOSE:
c Calculation of the spectrum of e-bremsstrahlung in the fully ionized,
c hydrogen-like, or helium-like gas as well as in neutral H, He, and 
c heavier in Born approximation; valid for E0,Egam >~ 0.01 MeV and
c beta0, beta1 >> 2*pi*Af*Z, corresponds to gamma=1.001
c REFERENCES:
c [1] Blumenthal & Gould 1970,Rev.Mod.Phys.42,237 (BG70); 
c [2] Gould 1969, Phys.Rev.185,72 (G69);
c [3] Koch & Motz 1959, Rev.Mod.Phys.31,920 (KM59).
c [4] Strong A.W., Moskalenko I.V., Reimer O. 2000, ApJ 537, 763
c INPUT/OUTPUT parameters:
c Egam, E0 [MeV] are the energy of a gamma and total electron energy,
c respectively;
c IZ1 is the nucleus charge; 
c Ne1 = 0 for unshielded charge, = 1 for 1-e atoms, = 2 for 2-e atoms.
c dSdK [barn/MeV] is the differential cross section (output).
c Note:
c #1 a one-parameter Hylleraas function is used for helium-like ions,
c    while Hartree-Fock functions are used for heutral He atoms;
c #2 a contribution of free electrons could be estimated by calling
c    the subroutine with Z=1, Ne=0 and multiplying dSdk by ionization level.
c***********************************************************************
    implicit real*8 (a-h,m,o-z)
*/
  /* this is in pp_meson.f :
      common/delta/ delta, IZ, Ne
     2      /Hartree/ DD(11),PH1(11),PH2(11)  ! DD = delta/(2.d0*Af)
c He formfactors
      data DD/0., 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10./,
     #     PH1/134.60, 133.85, 133.11, 130.86, 127.17, 
     #         120.35, 104.60,  89.94,  74.19,  54.26,  40.94/,
     #     PH2/131.40, 130.51, 130.33, 129.26, 126.76,
     #         120.80, 105.21,  89.46,  73.03,  51.84,  37.24/ 

      external FPHI1, FPHI2
  */
//////////////////////////////////////////////////////////////////////////////////////////////
