using namespace std;
#include <vector>
#include <valarray>
#include <string>
#include <iostream>
#include "CR_Spectrum.h"


void CR_Spectrum::init(string name_,valarray<double> parameters_, int debug_)
{
  cout<<">>CR_Spectrum::init"<<endl;
  name=name_;
  parameters.resize(parameters_.size());
  parameters=parameters_;

  debug=debug_;
  initialized=1;

  cout<<" name of CR spectrum :" <<name;

  cout<<"  parameters: ";
  for(int i=0;i<parameters.size();i++) cout<<parameters[i]<<" " ; cout<<endl;

  cout<<"<<CR_Spectrum::init"<<endl;
  return;
}

 /////////////////////////////////////////////////////////////////

// general CR spectra

valarray<double> CR_Spectrum::flux(valarray<double> E)
{
  if(debug==1)cout<<">>CR_Spectrum::flux name="<<name<<endl;

  valarray<double> flux_(E.size());
//////////////////////////////////////////////////////////////////

  if(name=="powerlaw")
  {
   double constant=parameters[1];
   double ge      =parameters[2];
   flux_          =constant*pow(E,-ge);
  }


//////////////////////////////////////////////////////////////////
  if(name=="Shikaze protons")
  {
    A=1.0;
    Z=1.0;

 // Shikaze et al 2007 http://adsabs.harvard.edu/abs/2007APh....28..154S
  // Fig 8 and 9 parameters
  // flux = A beta^P1 R^-P2

   double  Phi  =parameters[0]; // in MV as standard
   double  Phi0 = 591.;// the value used by Shikaze etal for 1998 data
   double dPhiGV=(Phi0-Phi)/1000.;// relative modulation parameter in GV 

  // protons
  double Ap=1.94e4*1e-4;// from m-2 sr-1 s-1 GeV-1 to cm-2  sr-1 s-1 GeV-1
  double P1=0.70;
  double P2=2.76;
  double m_p=0.938;     // proton mass in GeV

 
  valarray<double> Ek,g,beta,R;
  Ek  .resize(E.size());
  g   .resize(E.size());
  beta.resize(E.size());
  R   .resize(E.size());

  Ek   = E-m_p;                        // kinetic energy per nucleon

  if(dPhiGV==0.0) // no modulation correction, use LIS as given
 {
  g    = E/m_p;                        // Lorentz factor of proton
  beta = sqrt(1.0-1.0/g/g);            // velocity/c     of proton
  R    = sqrt(pow(E,2.)-pow(m_p,2.));  // rigidity = momentum

  flux_= Ap * pow(beta,P1) * pow(R,-P2);

  for(int i=0;i<E.size();i++)if(E[i]<m_p)flux_[i]=0.;

  for(int i=0;i<E.size();i++)
    {
      

     cout<<"E = "<< E[i]<< " Ek="<<Ek[i]<<" beta="<<beta[i]<<" R="<<R[i] 
                      <<" CR_Spectrum Shikaze proton spectrum="<< flux_[i]<<"  cm-2  sr-1 s-1 GeV-1"<<endl;

    }

   }//if

 else      // apply modulation correction to LIS
 {

/*
Gleeson & Axford,  ApJ 154, 1011 (1968)
The formulae here are based on equation (16) of that paper.

 J(E)           J(E + |Z| e phi)
________   =    _________________________ 
E^2-Eo^2         (E + |Z| e phi)^2 - Eo^2        

where E = total energy of particle, Eo = rest mass energy of particle.

see also galplot: modulate.cc
 */

  valarray<double>  E_prime, g_prime, beta_prime,R_prime,factor;
  E_prime    .resize(E.size());
  g_prime    .resize(E.size());
  beta_prime .resize(E.size());
  R_prime    .resize(E.size());
  factor     .resize(E.size());

  E_prime = E + Z *dPhiGV/A;                         // total energy per nucleon prime
  g_prime = E_prime/m_p;                             // Lorentz factor 
  beta_prime = sqrt(1.0-1.0/g_prime/g_prime);        // velocity/c    
  R_prime = sqrt(pow(A*E_prime,2.)-pow(A*m_p,2.))/Z; // rigidity = particle momentum/Z

  factor  = ( pow(A*E,2.) - pow(A*m_p,2.)) / ( pow(A*E_prime,2.) - pow(A*m_p,2.))   ;

  flux_   = Ap *pow(beta_prime,P1)*pow(R_prime,-P2)* factor;

  for(int i=0;i<E.size();i++)if(E_prime[i]<m_p)flux_[i]=0.;

  if(debug==1)
  for(int i=0;i<E.size();i++)
    cout<<"Phi="<<Phi<<" MV"<<" dPhiGV="<<dPhiGV<<" find E  = "<< E[i] <<" Ek="<<Ek[i]<<" GeV/n  E_prime="<<E_prime[i]
        <<" R_prime="<<R_prime[i]
        <<" CR_Spectrum Shikaze modulated proton spectrum="<< flux_[i]<<"  cm-2  sr-1 s-1 GeV-1" <<endl;


 }//else

  }

//////////////////////////////////////////////////////////////////

  if(name=="Shikaze Helium")
  {
    A=4.0;
    Z=2.0;

 // Shikaze et al 2007 http://adsabs.harvard.edu/abs/2007APh....28..154S
  // Fig 8 and 9 parameters
  // flux = A beta^P1 R^-P2

   double  Phi  =parameters[0]; // in MV as standard
   double  Phi0 = 591.;// the value used by Shikaze etal for 1998 data
   double dPhiGV=(Phi0-Phi)/1000.;// relative modulation parameter in GV 
 
  // Helium

   double     AHe=7.10e3*1e-4;// from m-2  sr-1 s-1 (GeV/n)-1 to cm-2  sr-1 s-1 (GeV/n)-1
   double      P1=0.50;
   double      P2=2.78;
   double     m_p=0.938;     // proton mass in GeV

  valarray<double> Ek,g,beta,R;
  Ek  .resize(E.size());
  g   .resize(E.size());
  beta.resize(E.size());
  R   .resize(E.size());

  Ek   = E-m_p;                              // kinetic energy per nucleon

 if(dPhiGV==0.0)// no modulation correction, use LIS as given
 {

  g    = E/m_p;                              // Lorentz factor
  beta = sqrt(1.0-1.0/g/g);                  // velocity/c    
  R    = sqrt(pow(A*E,2.)-pow(A*m_p,2.))/Z;  // rigidity = particle momentum/Z

  flux_= AHe * pow(beta,P1) * pow(R,-P2);

  for(int i=0;i<E.size();i++)if(E[i]<m_p)flux_[i]=0.;

  if(debug==1)
  for(int i=0;i<E.size();i++)
    {
      

     cout<<" find E = "<< E[i]<< " Ek="<<Ek[i]<<" beta="<<beta[i]<<" R="<<R[i] 
                      <<" CR_Spectrum Shikaze Helium spectrum="<< flux_[i]<<"  cm-2  sr-1 s-1 GeV-1"<<endl;

    }

  }//if

 else      // apply modulation correction to LIS
 {
 /*
Gleeson & Axford,  ApJ 154, 1011 (1968)
The formulae here are based on equation (16) of that paper.

 J(E)           J(E + |Z| e phi)
________   =    _________________________ 
E^2-Eo^2         (E + |Z| e phi)^2 - Eo^2        

where E = total energy of particle, Eo = rest mass energy of particle.

see also galplot: modulate.cc
 */
  valarray<double>  E_prime, g_prime, beta_prime,R_prime,factor;
  E_prime    .resize(E.size());
  g_prime    .resize(E.size());
  beta_prime .resize(E.size());
  R_prime    .resize(E.size());
  factor     .resize(E.size());

  E_prime = E + Z *dPhiGV/A;                         // total energy per nucleon prime
  g_prime = E_prime/m_p;                             // Lorentz factor 
  beta_prime = sqrt(1.0-1.0/g_prime/g_prime);        // velocity/c    
  R_prime = sqrt(pow(A*E_prime,2.)-pow(A*m_p,2.))/Z; // rigidity = particle momentum/Z

  factor  = ( pow(A*E,2.) - pow(A*m_p,2.)) / ( pow(A*E_prime,2.) - pow(A*m_p,2.))   ;

  flux_   = AHe*pow(beta_prime,P1)*pow(R_prime,-P2)* factor;

  if(debug==1)
  for(int i=0;i<E.size();i++)
    cout<<"Phi="<<Phi<<" MV"<<" dPhiGV="<<dPhiGV<<" find E  = "<< E[i] <<" Ek="<<Ek[i]<<" GeV/n  E_prime="<<E_prime[i]
        <<" R_prime="<<R_prime[i]
        <<" CR_Spectrum Shikaze modulated helium spectrum="<< flux_[i]<<"  cm-2  sr-1 s-1 GeV-1" <<endl;

  for(int i=0;i<E.size();i++)if(E_prime[i]<m_p)flux_[i]=0.;
 }//else



  }
  ///////////// Shikaze MeV versions


//////////////////////////////////////////////////////////////////
  if(name=="Shikaze protons MeV")
  {
    // same as Shikaze protons but interface all units in MeV: input E and output flux MeV-1
    A=1.0;
    Z=1.0;

 // Shikaze et al 2007 http://adsabs.harvard.edu/abs/2007APh....28..154S
  // Fig 8 and 9 parameters
  // flux = A beta^P1 R^-P2

   double  Phi  =parameters[0]; // in MV as standard
   double  Phi0 = 591.;// the value used by Shikaze etal for 1998 data
   double dPhiGV=(Phi0-Phi)/1000.;// relative modulation parameter in GV 

  // protons
  double Ap=1.94e4*1e-4;// from m-2 sr-1 s-1 GeV-1 to cm-2  sr-1 s-1 GeV-1
  double P1=0.70;
  double P2=2.76;
  double m_p=0.938;     // proton mass in GeV

  
 
  valarray<double> Ek,g,beta,R;
  Ek  .resize(E.size());
  g   .resize(E.size());
  beta.resize(E.size());
  R   .resize(E.size());

  valarray<double> EGeV;// since formulae in GeV
  EGeV  .resize(E.size());
  EGeV = E*1.e-3;// MeV->GeV

  


//Ek   = E   -m_p;                        // kinetic energy per nucleon
  Ek   = EGeV-m_p;                        // kinetic energy per nucleon

  if(dPhiGV==0.0) // no modulation correction, use LIS as given
 {
//g    = E   /m_p;                        // Lorentz factor of proton
  g    = EGeV/m_p;                        // Lorentz factor of proton
  beta = sqrt(1.0-1.0/g/g);            // velocity/c     of proton
//R    = sqrt(pow(E,   2.)-pow(m_p,2.));  // rigidity = momentum
  R    = sqrt(pow(EGeV,2.)-pow(m_p,2.));  // rigidity = momentum

  flux_= Ap * pow(beta,P1) * pow(R,-P2);

  flux_ *= 1.e-3;// GeV-1 -> MeV-1

  for(int i=0;i<E.size();i++)if(EGeV[i]<m_p)flux_[i]=0.;

  if(debug==1)
  for(int i=0;i<E.size();i++)
    {
      

     cout<<"E = "<< E[i]<< " Ek="<<Ek[i]<<" beta="<<beta[i]<<" R="<<R[i] 
                      <<" CR_Spectrum Shikaze proton spectrum MeV="<< flux_[i]<<"  cm-2  sr-1 s-1 MeV-1"<<endl;

    }

   }//if

 else      // apply modulation correction to LIS
 {

/*
Gleeson & Axford,  ApJ 154, 1011 (1968)
The formulae here are based on equation (16) of that paper.

 J(E)           J(E + |Z| e phi)
________   =    _________________________ 
E^2-Eo^2         (E + |Z| e phi)^2 - Eo^2        

where E = total energy of particle, Eo = rest mass energy of particle.

see also galplot: modulate.cc
 */

  valarray<double>  E_prime, g_prime, beta_prime,R_prime,factor;
  E_prime    .resize(E.size());
  g_prime    .resize(E.size());
  beta_prime .resize(E.size());
  R_prime    .resize(E.size());
  factor     .resize(E.size());

//E_prime = E    + Z *dPhiGV/A;                         // total energy per nucleon prime
  E_prime = EGeV + Z *dPhiGV/A;                         // total energy per nucleon prime
  g_prime = E_prime/m_p;                             // Lorentz factor 
  beta_prime = sqrt(1.0-1.0/g_prime/g_prime);        // velocity/c    
  R_prime = sqrt(pow(A*E_prime,2.)-pow(A*m_p,2.))/Z; // rigidity = particle momentum/Z

//factor  = ( pow(A*E   ,2.) - pow(A*m_p,2.)) / ( pow(A*E_prime,2.) - pow(A*m_p,2.))   ;
  factor  = ( pow(A*EGeV,2.) - pow(A*m_p,2.)) / ( pow(A*E_prime,2.) - pow(A*m_p,2.))   ;

  flux_   = Ap *pow(beta_prime,P1)*pow(R_prime,-P2)* factor;

  flux_ *= 1.e-3;// GeV-1 -> MeV-1

  for(int i=0;i<E.size();i++)if(E_prime[i]<m_p)flux_[i]=0.;

  if(debug==1)
  for(int i=0;i<E.size();i++)
    cout<<"Phi="<<Phi<<" MV"<<" dPhiGV="<<dPhiGV<<" find E  = "<< E[i]<<" MeV/n" 
        <<" Ek="<<Ek[i]<<" GeV/n  E_prime="<<E_prime[i]
        <<" R_prime="<<R_prime[i]
        <<" CR_Spectrum Shikaze modulated proton spectrum MeV="<< flux_[i]<<"  cm-2  sr-1 s-1 MeV-1" <<endl;


 }//else



  }

//////////////////////////////////////////////////////////////////

  if(name=="Shikaze Helium MeV")
  {
   // same as Shikaze Helium  but interface all units in MeV: input E and output flux MeV-1
    A=4.0;
    Z=2.0;

 // Shikaze et al 2007 http://adsabs.harvard.edu/abs/2007APh....28..154S
  // Fig 8 and 9 parameters
  // flux = A beta^P1 R^-P2

   double  Phi  =parameters[0]; // in MV as standard
   double  Phi0 = 591.;// the value used by Shikaze etal for 1998 data
   double dPhiGV=(Phi0-Phi)/1000.;// relative modulation parameter in GV 
 
  // Helium

   double     AHe=7.10e3*1e-4;// from m-2  sr-1 s-1 (GeV/n)-1 to cm-2  sr-1 s-1 (GeV/n)-1
   double      P1=0.50;
   double      P2=2.78;
   double     m_p=0.938;     // proton mass in GeV

  valarray<double> Ek,g,beta,R;
  Ek  .resize(E.size());
  g   .resize(E.size());
  beta.resize(E.size());
  R   .resize(E.size());

  valarray<double> EGeV;// since formulae in GeV
  EGeV  .resize(E.size());
  EGeV = E*1.e-3;// MeV->GeV

//Ek   = E   -m_p;                              // kinetic energy per nucleon
  Ek   = EGeV-m_p;                              // kinetic energy per nucleon

 if(dPhiGV==0.0)// no modulation correction, use LIS as given
 {

//g    = E   /m_p;                              // Lorentz factor
  g    = EGeV/m_p;                              // Lorentz factor
  beta = sqrt(1.0-1.0/g/g);                  // velocity/c    
//R    = sqrt(pow(A*E   ,2.)-pow(A*m_p,2.))/Z;  // rigidity = particle momentum/Z
  R    = sqrt(pow(A*EGeV,2.)-pow(A*m_p,2.))/Z;  // rigidity = particle momentum/Z

  flux_= AHe * pow(beta,P1) * pow(R,-P2);
  flux_ *= 1.e-3;// GeV-1 -> MeV-1

  for(int i=0;i<E.size();i++)if(EGeV[i]<m_p)flux_[i]=0.;

  if(debug==1)
  for(int i=0;i<E.size();i++)
    {
      

     cout<<" find E = "<< E[i]<< "MeV Ek="<<Ek[i]<<"GeV  beta="<<beta[i]<<" R="<<R[i] 
                      <<" CR_Spectrum Shikaze Helium spectrum MeV="<< flux_[i]<<"  cm-2  sr-1 s-1 MeV-1"<<endl;

    }

  }//if

 else      // apply modulation correction to LIS
 {
 /*
Gleeson & Axford,  ApJ 154, 1011 (1968)
The formulae here are based on equation (16) of that paper.

 J(E)           J(E + |Z| e phi)
________   =    _________________________ 
E^2-Eo^2         (E + |Z| e phi)^2 - Eo^2        

where E = total energy of particle, Eo = rest mass energy of particle.

see also galplot: modulate.cc
 */
  valarray<double>  E_prime, g_prime, beta_prime,R_prime,factor;
  E_prime    .resize(E.size());
  g_prime    .resize(E.size());
  beta_prime .resize(E.size());
  R_prime    .resize(E.size());
  factor     .resize(E.size());

//E_prime = E    + Z *dPhiGV/A;                         // total energy per nucleon prime
  E_prime = EGeV + Z *dPhiGV/A;                         // total energy per nucleon prime
  g_prime = E_prime/m_p;                             // Lorentz factor 
  beta_prime = sqrt(1.0-1.0/g_prime/g_prime);        // velocity/c    
  R_prime = sqrt(pow(A*E_prime,2.)-pow(A*m_p,2.))/Z; // rigidity = particle momentum/Z

//factor  = ( pow(A*E,   2.) - pow(A*m_p,2.)) / ( pow(A*E_prime,2.) - pow(A*m_p,2.))   ;
  factor  = ( pow(A*EGeV,2.) - pow(A*m_p,2.)) / ( pow(A*E_prime,2.) - pow(A*m_p,2.))   ;

  flux_   = AHe*pow(beta_prime,P1)*pow(R_prime,-P2)* factor;
  flux_ *= 1.e-3; // GeV-1 -> MeV-1

  if(debug==1)
  for(int i=0;i<E.size();i++)
    cout<<"Phi="<<Phi<<" MV"<<" dPhiGV="<<dPhiGV<<" find E  = "<< E[i] <<"MeV Ek="<<Ek[i]<<" GeV/n  E_prime="<<E_prime[i]
        <<" R_prime="<<R_prime[i]
        <<" CR_Spectrum Shikaze modulated helium spectrum MeV="<< flux_[i]<<"  cm-2  sr-1 s-1 MeV-1" <<endl;

  for(int i=0;i<E.size();i++)if(E_prime[i]<m_p)flux_[i]=0.;

 }//else



  }

  ////// Shikaze MeV kinetic energy

//////////////////////////////////////////////////////////////////
  if(name=="Shikaze protons kinetic energy MeV")
  {
    // same as Shikaze protons but interface all units in MeV: input E and output flux MeV-1
    A=1.0;
    Z=1.0;

 // Shikaze et al 2007 http://adsabs.harvard.edu/abs/2007APh....28..154S
  // Fig 8 and 9 parameters
  // flux = A beta^P1 R^-P2

   double  Phi  =parameters[0]; // in MV as standard
   double  Phi0 = 591.;// the value used by Shikaze etal for 1998 data
   double dPhiGV=(Phi0-Phi)/1000.;// relative modulation parameter in GV 

  // protons
  double Ap=1.94e4*1e-4;// from m-2 sr-1 s-1 GeV-1 to cm-2  sr-1 s-1 GeV-1
  double P1=0.70;
  double P2=2.76;
  double m_p=0.938;     // proton mass in GeV

 
  valarray<double> Ek,g,beta,R;
  Ek  .resize(E.size());
  g   .resize(E.size());
  beta.resize(E.size());
  R   .resize(E.size());

  valarray<double> EGeV;// since formulae in GeV
  EGeV  .resize(E.size());
  EGeV = E*1.e-3;// MeV->GeV

//Ek   = E   -m_p;                        // kinetic energy per nucleon

///  Ek   = EGeV-m_p;                        // kinetic energy per nucleon: as input in this case

    Ek   = EGeV;
    EGeV += m_p; // total energy from kinetic energy

  if(dPhiGV==0.0) // no modulation correction, use LIS as given
 {
//g    = E   /m_p;                        // Lorentz factor of proton
  g    = EGeV/m_p;                        // Lorentz factor of proton
  beta = sqrt(1.0-1.0/g/g);            // velocity/c     of proton
//R    = sqrt(pow(E,   2.)-pow(m_p,2.));  // rigidity = momentum
  R    = sqrt(pow(EGeV,2.)-pow(m_p,2.));  // rigidity = momentum

  flux_= Ap * pow(beta,P1) * pow(R,-P2);

  flux_ *= 1.e-3;// GeV-1 -> MeV-1

  if(debug==1)
  for(int i=0;i<E.size();i++)
    {
      

     cout<<"E = "<< E[i]<< " Ek="<<Ek[i]<<" beta="<<beta[i]<<" R="<<R[i] 
                      <<" CR_Spectrum Shikaze proton spectrum MeV="<< flux_[i]<<"  cm-2  sr-1 s-1 MeV-1"<<endl;

    }

   }//if

 else      // apply modulation correction to LIS
 {

/*
Gleeson & Axford,  ApJ 154, 1011 (1968)
The formulae here are based on equation (16) of that paper.

 J(E)           J(E + |Z| e phi)
________   =    _________________________ 
E^2-Eo^2         (E + |Z| e phi)^2 - Eo^2        

where E = total energy of particle, Eo = rest mass energy of particle.

see also galplot: modulate.cc
 */

  valarray<double>  E_prime, g_prime, beta_prime,R_prime,factor;
  E_prime    .resize(E.size());
  g_prime    .resize(E.size());
  beta_prime .resize(E.size());
  R_prime    .resize(E.size());
  factor     .resize(E.size());

//E_prime = E    + Z *dPhiGV/A;                         // total energy per nucleon prime
  E_prime = EGeV + Z *dPhiGV/A;                         // total energy per nucleon prime
  g_prime = E_prime/m_p;                             // Lorentz factor 
  beta_prime = sqrt(1.0-1.0/g_prime/g_prime);        // velocity/c    
  R_prime = sqrt(pow(A*E_prime,2.)-pow(A*m_p,2.))/Z; // rigidity = particle momentum/Z

//factor  = ( pow(A*E   ,2.) - pow(A*m_p,2.)) / ( pow(A*E_prime,2.) - pow(A*m_p,2.))   ;
  factor  = ( pow(A*EGeV,2.) - pow(A*m_p,2.)) / ( pow(A*E_prime,2.) - pow(A*m_p,2.))   ;

  flux_   = Ap *pow(beta_prime,P1)*pow(R_prime,-P2)* factor;

  flux_ *= 1.e-3;// GeV-1 -> MeV-1

  if(debug==1)
  for(int i=0;i<E.size();i++)
    cout<<"Phi="<<Phi<<" MV"<<" dPhiGV="<<dPhiGV<<" find E  = "<< E[i]<<" MeV/n" 
        <<" Ek="<<Ek[i]<<" GeV/n  E_prime="<<E_prime[i]
        <<" R_prime="<<R_prime[i]
        <<" CR_Spectrum Shikaze modulated proton spectrum MeV="<< flux_[i]<<"  cm-2  sr-1 s-1 MeV-1" <<endl;


 }//else

  }

//////////////////////////////////////////////////////////////////

  if(name=="Shikaze Helium kinetic energy MeV")
  {
   // same as Shikaze Helium  but interface all units in MeV: input E and output flux MeV-1
    A=4.0;
    Z=2.0;

 // Shikaze et al 2007 http://adsabs.harvard.edu/abs/2007APh....28..154S
  // Fig 8 and 9 parameters
  // flux = A beta^P1 R^-P2

   double  Phi  =parameters[0]; // in MV as standard
   double  Phi0 = 591.;// the value used by Shikaze etal for 1998 data
   double dPhiGV=(Phi0-Phi)/1000.;// relative modulation parameter in GV 
 
  // Helium

   double     AHe=7.10e3*1e-4;// from m-2  sr-1 s-1 (GeV/n)-1 to cm-2  sr-1 s-1 (GeV/n)-1
   double      P1=0.50;
   double      P2=2.78;
   double     m_p=0.938;     // proton mass in GeV

  valarray<double> Ek,g,beta,R;
  Ek  .resize(E.size());
  g   .resize(E.size());
  beta.resize(E.size());
  R   .resize(E.size());

  valarray<double> EGeV;// since formulae in GeV
  EGeV  .resize(E.size());
  EGeV = E*1.e-3;// MeV->GeV

//Ek   = E   -m_p;                              // kinetic energy per nucleon
///  Ek   = EGeV-m_p;                              // kinetic energy per nucleon
    Ek = EGeV;// input in kinetic energy per nucleon in this case

    EGeV += m_p; // total energy/nucleon from kinetic energy/nucleon

 if(dPhiGV==0.0)// no modulation correction, use LIS as given
 {

//g    = E   /m_p;                              // Lorentz factor
  g    = EGeV/m_p;                              // Lorentz factor
  beta = sqrt(1.0-1.0/g/g);                  // velocity/c    
//R    = sqrt(pow(A*E   ,2.)-pow(A*m_p,2.))/Z;  // rigidity = particle momentum/Z
  R    = sqrt(pow(A*EGeV,2.)-pow(A*m_p,2.))/Z;  // rigidity = particle momentum/Z

  flux_= AHe * pow(beta,P1) * pow(R,-P2);
  flux_ *= 1.e-3;// GeV-1 -> MeV-1

  if(debug==1)
  for(int i=0;i<E.size();i++)
    {
      

     cout<<" find E = "<< E[i]<< "MeV Ek="<<Ek[i]<<"GeV  beta="<<beta[i]<<" R="<<R[i] 
                      <<" CR_Spectrum Shikaze Helium spectrum MeV="<< flux_[i]<<"  cm-2  sr-1 s-1 MeV-1"<<endl;

    }

  }//if

 else      // apply modulation correction to LIS
 {
 /*
Gleeson & Axford,  ApJ 154, 1011 (1968)
The formulae here are based on equation (16) of that paper.

 J(E)           J(E + |Z| e phi)
________   =    _________________________ 
E^2-Eo^2         (E + |Z| e phi)^2 - Eo^2        

where E = total energy of particle, Eo = rest mass energy of particle.

see also galplot: modulate.cc
 */
  valarray<double>  E_prime, g_prime, beta_prime,R_prime,factor;
  E_prime    .resize(E.size());
  g_prime    .resize(E.size());
  beta_prime .resize(E.size());
  R_prime    .resize(E.size());
  factor     .resize(E.size());

//E_prime = E    + Z *dPhiGV/A;                         // total energy per nucleon prime
  E_prime = EGeV + Z *dPhiGV/A;                         // total energy per nucleon prime
  g_prime = E_prime/m_p;                             // Lorentz factor 
  beta_prime = sqrt(1.0-1.0/g_prime/g_prime);        // velocity/c    
  R_prime = sqrt(pow(A*E_prime,2.)-pow(A*m_p,2.))/Z; // rigidity = particle momentum/Z

//factor  = ( pow(A*E,   2.) - pow(A*m_p,2.)) / ( pow(A*E_prime,2.) - pow(A*m_p,2.))   ;
  factor  = ( pow(A*EGeV,2.) - pow(A*m_p,2.)) / ( pow(A*E_prime,2.) - pow(A*m_p,2.))   ;

  flux_   = AHe*pow(beta_prime,P1)*pow(R_prime,-P2)* factor;
  flux_ *= 1.e-3; // GeV-1 -> MeV-1

  if(debug==1)
  for(int i=0;i<E.size();i++)
    cout<<"Phi="<<Phi<<" MV"<<" dPhiGV="<<dPhiGV<<" find E  = "<< E[i] <<"MeV Ek="<<Ek[i]<<" GeV/n  E_prime="<<E_prime[i]
        <<" R_prime="<<R_prime[i]
        <<" CR_Spectrum Shikaze modulated helium spectrum MeV="<< flux_[i]<<"  cm-2  sr-1 s-1 MeV-1" <<endl;


 }//else

  }
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
  if(name=="Gaisser protons")
  {
    A=1.0;
    Z=1.0;

 //  Gaisser etal 2001 fits to heliospheric spectra away from solar max
 // http://adsabs.harvard.edu/abs/2001ICRC....5.1643G
 // Table 1 parameters
 // flux = 

   double  Phi  =parameters[0]; // in MV as standard. The modulation which will be assumed to demodulate the spectrum.
   double  Phi0 = 0.; // start heliospheric
   double dPhiGV=(Phi0-Phi)/1000.;// relative modulation parameter in GV : here negative since demodulation

  // protons
   

  double alpha = 2.74;
  double K     = 14900.;//  m-2  sr-1 s-1 GeV-1
  double     b = 2.15;
  double     c = 0.21;

  double     m_p=0.938;     // proton mass in GeV

  valarray<double> Ek;
  Ek  .resize(E.size());
 
 
  Ek   = E-m_p;                        // kinetic energy per nucleon

  if(dPhiGV==0.0) // no modulation correction, use LIS as given
 {
 

   flux_= K*1.0e-4 * pow(Ek + b*exp(-c*sqrt(Ek)),-alpha); // equation (1) ->  cm-2  sr-1 s-1 GeV-1

  for(int i=0;i<E.size();i++)
    {
      

     cout<<"E = "<< E[i]<< " Ek="<<Ek[i]
                      <<" CR_Spectrum Gaisser proton spectrum="<< flux_[i]<<"  cm-2  sr-1 s-1 GeV-1"<<endl;

    }

   }//if

 else      // apply demodulation correction 
 {

/*
Gleeson & Axford,  ApJ 154, 1011 (1968)
The formulae here are based on equation (16) of that paper.

 J(E)           J(E + |Z| e phi)
________   =    _________________________ 
E^2-Eo^2         (E + |Z| e phi)^2 - Eo^2        

where E = total energy of particle, Eo = rest mass energy of particle.

see also galplot: modulate.cc
 */

  valarray<double>  E_prime,Ek_prime, factor;

  E_prime    .resize(E.size());
  Ek_prime   .resize(E.size());
  factor     .resize(E.size());

  Ek_prime = Ek + Z *dPhiGV/A;                         // kinetic energy per nucleon prime
  E_prime  = E  + Z *dPhiGV/A;                         // total energy per nucleon prime
  

  factor  = ( pow(A*E,2.) - pow(A*m_p,2.)) / ( pow(A*E_prime,2.) - pow(A*m_p,2.))   ;

  flux_= K*1.0e-4 * pow(Ek_prime + b*exp(-c*sqrt(Ek_prime)),-alpha); // equation (1) ->  cm-2  sr-1 s-1 GeV-1* factor;
  flux_*=factor;


  for(int i=0;i<E.size();i++)if(Ek_prime[i]<=0.)flux_[i]=0.; // to avoid nans

  for(int i=0;i<E.size();i++)
    cout<<"Phi="<<Phi<<" MV"<<" dPhiGV="<<dPhiGV<<" E total  = "<< E[i] <<" Ek="<<Ek[i]<<" GeV/n  Ek_prime="<<Ek_prime[i]
        <<" GeV/n  E_prime="<<E_prime[i]
        <<" CR_Spectrum Gaisser demodulated proton spectrum="<< flux_[i]<<"  cm-2  sr-1 s-1 GeV-1" <<endl;


 }//else

  }// Gaisser protons



//////////////////////////////////////////////////////////////////
  if(name=="Gaisser Helium high")
  {
    A=4.0;
    Z=2.0;

 //  Gaisser etal 2001 fits to heliospheric spectra away from solar max
 // http://adsabs.harvard.edu/abs/2001ICRC....5.1643G
 // Table 1 parameters
 

   double  Phi  =parameters[0]; // in MV as standard. The modulation which will be assumed to demodulate the spectrum.
   double  Phi0 = 0.; // start heliospheric
   double dPhiGV=(Phi0-Phi)/1000.;// relative modulation parameter in GV : here negative since demodulation

  // Helium high
   

  double alpha = 2.64;
  double K     =   600.;//  m-2  sr-1 s-1 GeV-1
  double     b = 1.25;
  double     c = 0.14;

  double     m_p=0.938;     // proton mass in GeV

  valarray<double> Ek;
  Ek  .resize(E.size());
 
 
  Ek   = E-m_p;                        // kinetic energy per nucleon

  if(dPhiGV==0.0) // no modulation correction, use LIS as given
 {
 

   flux_= K*1.0e-4 * pow(Ek + b*exp(-c*sqrt(Ek)),-alpha); // equation (1) ->  cm-2  sr-1 s-1 GeV-1

  for(int i=0;i<E.size();i++)
    {
      

     cout<<"E = "<< E[i]<< " Ek="<<Ek[i]
                      <<" CR_Spectrum Gaisser Helium high spectrum="<< flux_[i]<<"  cm-2  sr-1 s-1 GeV-1"<<endl;

    }

   }//if

 else      // apply demodulation correction 
 {

/*
Gleeson & Axford,  ApJ 154, 1011 (1968)
The formulae here are based on equation (16) of that paper.

 J(E)           J(E + |Z| e phi)
________   =    _________________________ 
E^2-Eo^2         (E + |Z| e phi)^2 - Eo^2        

where E = total energy of particle, Eo = rest mass energy of particle.

see also galplot: modulate.cc
 */

  valarray<double>  E_prime,Ek_prime, factor;

  E_prime    .resize(E.size());
  Ek_prime   .resize(E.size());
  factor     .resize(E.size());

  Ek_prime = Ek + Z *dPhiGV/A;                         // kinetic energy per nucleon prime
  E_prime  = E  + Z *dPhiGV/A;                         // total energy per nucleon prime
  

  factor  = ( pow(A*E,2.) - pow(A*m_p,2.)) / ( pow(A*E_prime,2.) - pow(A*m_p,2.))   ;

  flux_= K*1.0e-4 * pow(Ek_prime + b*exp(-c*sqrt(Ek_prime)),-alpha); // equation (1) ->  cm-2  sr-1 s-1 GeV-1* factor;
  flux_*=factor;


  for(int i=0;i<E.size();i++)if(Ek_prime[i]<=0.)flux_[i]=0.; // to avoid nans

  for(int i=0;i<E.size();i++)
    cout<<"Phi="<<Phi<<" MV"<<" dPhiGV="<<dPhiGV<<" E total  = "<< E[i] <<" Ek="<<Ek[i]<<" GeV/n  Ek_prime="<<Ek_prime[i]
        <<" GeV/n  E_prime="<<E_prime[i]
        <<" CR_Spectrum Gaisser demodulated Helium high spectrum="<< flux_[i]<<"  cm-2  sr-1 s-1 GeV-1" <<endl;


 }//else

  }// Gaisser Helium high

//////////////////////////////////////////////////////////////////
  if(name=="Gaisser Helium low")
  {
    A=4.0;
    Z=2.0;

 //  Gaisser etal 2001 fits to heliospheric spectra away from solar max
 // http://adsabs.harvard.edu/abs/2001ICRC....5.1643G
 // Table 1 parameters
 

   double  Phi  =parameters[0]; // in MV as standard. The modulation which will be assumed to demodulate the spectrum.
   double  Phi0 = 0.; // start heliospheric
   double dPhiGV=(Phi0-Phi)/1000.;// relative modulation parameter in GV : here negative since demodulation

  // Helium low
   

  double alpha = 2.74;
  double K     =   750.;//  m-2  sr-1 s-1 GeV-1
  double     b = 1.50;
  double     c = 0.30;

  double     m_p=0.938;     // proton mass in GeV

  valarray<double> Ek;
  Ek  .resize(E.size());
 
 
  Ek   = E-m_p;                        // kinetic energy per nucleon

  if(dPhiGV==0.0) // no modulation correction, use LIS as given
 {
 

   flux_= K*1.0e-4 * pow(Ek + b*exp(-c*sqrt(Ek)),-alpha); // equation (1) ->  cm-2  sr-1 s-1 GeV-1

  for(int i=0;i<E.size();i++)
    {
      

     cout<<"E = "<< E[i]<< " Ek="<<Ek[i]
                      <<" CR_Spectrum Gaisser Helium spectrum low="<< flux_[i]<<"  cm-2  sr-1 s-1 GeV-1"<<endl;

    }

   }//if

 else      // apply demodulation correction 
 {

/*
Gleeson & Axford,  ApJ 154, 1011 (1968)
The formulae here are based on equation (16) of that paper.

 J(E)           J(E + |Z| e phi)
________   =    _________________________ 
E^2-Eo^2         (E + |Z| e phi)^2 - Eo^2        

where E = total energy of particle, Eo = rest mass energy of particle.

see also galplot: modulate.cc
 */

  valarray<double>  E_prime,Ek_prime, factor;

  E_prime    .resize(E.size());
  Ek_prime   .resize(E.size());
  factor     .resize(E.size());

  Ek_prime = Ek + Z *dPhiGV/A;                         // kinetic energy per nucleon prime
  E_prime  = E  + Z *dPhiGV/A;                         // total energy per nucleon prime
  

  factor  = ( pow(A*E,2.) - pow(A*m_p,2.)) / ( pow(A*E_prime,2.) - pow(A*m_p,2.))   ;

  flux_= K*1.0e-4 * pow(Ek_prime + b*exp(-c*sqrt(Ek_prime)),-alpha); // equation (1) ->  cm-2  sr-1 s-1 GeV-1* factor;
  flux_*=factor;


  for(int i=0;i<E.size();i++)if(Ek_prime[i]<=0.)flux_[i]=0.; // to avoid nans

  for(int i=0;i<E.size();i++)
    cout<<"Phi="<<Phi<<" MV"<<" dPhiGV="<<dPhiGV<<" E total  = "<< E[i] <<" Ek="<<Ek[i]<<" GeV/n  Ek_prime="<<Ek_prime[i]
        <<" GeV/n  E_prime="<<E_prime[i]
        <<" CR_Spectrum Gaisser demodulated Helium low spectrum="<< flux_[i]<<"  cm-2  sr-1 s-1 GeV-1" <<endl;


 }//else

  }// Gaisser Helium low


  /////////////////////////// final part

  if(debug==1)
    for(int i=0;i<E.size();i++) cout<<"CR_Spectrum: name="<<name<<" Z="<<Z<<" A="<<A
                                    <<"  E= "<<E[i]<<" flux="<<flux_[i]<<endl;

  if(debug==1)  cout<<"<<CR_Spectrum::flux"<<endl;

  return flux_;
}

////////////////////////////////////////////////////////////////////////////
// electron and positron spectra as function of distance from sun
////////////////////////////////////////////////////////////////////////////


///////////////////////////////
valarray<double> CR_Spectrum::fluxe(valarray<double> Ee, double r)
{
  if(debug==1)cout<<">>CR_Spectrum::fluxe name="<<name<<endl;

  valarray<double> result(Ee.size());
//////////////////////////////////////////////////////////////////
  if(name=="powerlaw")
  {
   double constant=parameters[1];
   double ge      =parameters[2];
   result=constant*pow(Ee,-ge);
  }
//////////////////////////////////////////////////////////////////
 if(name=="Fermi2011_electrons")
 { // Ackermann etal Fermi-LAT arXiv:1109.052 electrons fit, converting from GeV m^2 to MeV cm^2 
  double ge=3.19;
  result = 2.07e-9* pow(Ee/2.e4, -ge);
 }
//////////////////////////////////////////////////////////////////
 if(name=="Fermi2011_positrons")
 { // Ackermann etal Fermi-LAT arXiv:1109.052 positrons fit, converting from GeV m^2 to MeV cm^2 
  double ge=2.77;
  result = 2.02e-10* pow(Ee/2.e4, -ge);
 }

 //////////////////////////////////////////////////////////////////
 if(name=="Fermi2011_electrons_positrons")
 { // Ackermann etal Fermi-LAT arXiv:1109.052 electrons and positrons fits summed, converting from GeV m^2 to MeV cm^2
 
   double Ee_break=parameters[1];// electron break energy
   double ge1     =parameters[2];// electron index below break
   double Ep_break=parameters[3];// positron break energy
   double gp1     =parameters[4];// positron index below break

   double electrons_factor=parameters[5]; // electrons factor relative to reference
   double positrons_factor=parameters[6]; // positrons factor relative to reference

   if (debug==1) cout<<"Ee_break="<<Ee_break<<" ge1="<<ge1<<" Ep_break="<<Ep_break<<" gp1="<<gp1<<endl;

   valarray<double> electrons(Ee.size());
   valarray<double> positrons(Ee.size());

  // electrons
  double ge=3.19;
  electrons  = 2.07e-9 * pow(Ee/2.e4, -ge);


  for(int i=0;i<Ee.size();i++)  if (Ee[i]<Ee_break)
   {            electrons[i] *= pow(Ee[i]/Ee_break, -(ge1-ge));
     if (debug==1)
      cout<<Ee[i]<<" electron energy below break "<<" factor="<<pow(Ee[i]/Ee_break, -(ge1-ge))
          <<" electrons="<<electrons[i]<<endl;
   }

  // positrons
  double gp=2.77;
  positrons= 2.02e-10* pow(Ee/2.e4, -gp);

  for(int i=0;i<Ee.size();i++)  if (Ee[i]<Ep_break)
   {            positrons[i] *= pow(Ee[i]/Ep_break, -(gp1-gp));
     if (debug==1)
      cout<<Ee[i]<<" positron energy below break "<<" factor="<<pow(Ee[i]/Ep_break, -(gp1-gp))
          <<" positrons="<<positrons[i]<<endl;
   }

  result = electrons * electrons_factor + positrons * positrons_factor;
 }

 /////////////////////////////////////////////////////////////////

if(name=="Fermi2011_electrons_positrons_modulated1")
 { // Ackermann etal Fermi-LAT arXiv:1109.052 electrons and positrons fits summed, converting from GeV m^2 to MeV cm^2 
   // the ubiquituous force-field approximation
   // expressions from Moskalenko et al. ApJ 652, L65 (2006), Abdo etal 2011 ApJ 734, 116 (2011) equations (2) and (3)
   // Phi1: cycles 20/22  Phi2:cycle 21
 
   
   double Phi0    =parameters[0];

   double Phi1,Phi2;
   double rb=100.;
   double r0=10.;

   if(r>=r0)Phi1=Phi0/1.88*          (pow(r,-0.4)-pow(rb,-0.4)); 
   if(r< r0)Phi1=Phi0/1.88*(0.24+8.0*(pow(r,-0.1)-pow(r0,-0.1)));

   double Phi=Phi1-Phi0; // since using spectrum at 1AU, modulation is relative to this

  if(debug==1)cout<<"CR_Spectrum name="<<name<<" r="<<r<<" Phi="<<Phi<<endl;




   double Ee_break=parameters[1];// electron break energy
   double ge1     =parameters[2];// electron index below break
   double Ep_break=parameters[3];// positron break energy
   double gp1     =parameters[4];// positron index below break

   double electrons_factor=parameters[5]; // electrons factor relative to reference
   double positrons_factor=parameters[6]; // positrons factor relative to reference

   if (debug==1) cout<<"Ee_break="<<Ee_break<<" ge1="<<ge1<<" Ep_break="<<Ep_break<<" gp1="<<gp1<<endl;

   valarray<double> electrons(Ee.size());
   valarray<double> positrons(Ee.size());

   valarray<double> Eeprime=Ee+Phi;

  // electrons
  double ge = 3.19;
  electrons = 0.;
  //  electrons  = 2.07e-9 * pow(Eeprime/2.e4, -ge);

  // electron mininum energy applied everywhere, required r>1AU where Eeprime<Ee.
  for(int i=0;i<Ee.size();i++)    if (Eeprime[i]>Ee[0]) electrons[i]  = 2.07e-9 * pow(Eeprime[i]/2.e4, -ge);

  for(int i=0;i<Ee.size();i++)    if (Eeprime[i] < Ee_break && Eeprime[i]>Ee[0])
    {            electrons[i] *= pow( Eeprime[i] / Ee_break, -(ge1-ge));

     if (debug==1)
       cout<<" Ee="<<Ee[i]<< " Phi="<<Phi<<" electron energy + Phi below break="<<Eeprime[i]
           <<" factor="<<pow(Eeprime[i]/Ee_break, -(ge1-ge))
           <<" electrons="<<electrons[i]<<endl;
   }

  // positrons
  double gp = 2.77;
  positrons = 0;
  //  positrons= 2.02e-10* pow(Eeprime /2.e4, -gp);

  for(int i=0;i<Ee.size();i++)    if (Eeprime[i]>Ee[0])   positrons[i]= 2.02e-10* pow(Eeprime[i] /2.e4, -gp);

  for(int i=0;i<Ee.size();i++)    if (Eeprime[i] < Ep_break && Eeprime[i]>Ee[0] )
    {            positrons[i] *= pow( Eeprime[i] / Ep_break, -(gp1-gp));

     if (debug==1)
       cout<<" Ee="<<Ee[i]<< " Phi="<<Phi<<" positron energy + Phi below break="<<Eeprime[i]
           <<" factor="<<pow(Eeprime[i]/Ep_break, -(gp1-gp))
           <<" positrons="<<positrons[i]<<endl;
   }

  result = electrons * electrons_factor + positrons * positrons_factor;



  result *= pow(Ee/(Ee+Phi),2.0);
 }

////////////////////////////////////////////////////////////////

if(name=="Fermi2011_electrons_positrons_modulated2")
 { // Ackermann etal Fermi-LAT arXiv:1109.052 electrons and positrons fits summed, converting from GeV m^2 to MeV cm^2 
   // the ubiquituous force-field approximation
   // expressions from Moskalenko et al. ApJ 652, L65 (2006), Abdo etal 2011 ApJ 734, 116 (2011) equations (2) and (3)
   // Phi1: cycles 20/22  Phi2:cycle 21
 
   double Phi0    =parameters[0];

   double Phi1,Phi2;
   double rb=100.;
   double r0=10.;



  Phi2=Phi0     *          (pow(r,-0.1)-pow(rb,-0.1))/(1.0-pow(rb,-0.1));

  double Phi=Phi2-Phi0; // since using spectrum at 1AU, modulation is relative to this

  if(debug==1)cout<<"CR_Spectrum name="<<name<<" r="<<r<<" Phi="<<Phi<<endl;




   double Ee_break=parameters[1];// electron break energy
   double ge1     =parameters[2];// electron index below break
   double Ep_break=parameters[3];// positron break energy
   double gp1     =parameters[4];// positron index below break

   double electrons_factor=parameters[5]; // electrons factor relative to reference
   double positrons_factor=parameters[6]; // positrons factor relative to reference

   if (debug==1) cout<<"Ee_break="<<Ee_break<<" ge1="<<ge1<<" Ep_break="<<Ep_break<<" gp1="<<gp1<<endl;

   valarray<double> electrons(Ee.size());
   valarray<double> positrons(Ee.size());

   valarray<double> Eeprime=Ee+Phi;

  // electrons
  double ge = 3.19;
  electrons = 0.;
  //  electrons  = 2.07e-9 * pow(Eeprime/2.e4, -ge);

  // electron mininum energy applied everywhere, required r>1AU where Eeprime<Ee.
  for(int i=0;i<Ee.size();i++)    if (Eeprime[i]>Ee[0]) electrons[i]  = 2.07e-9 * pow(Eeprime[i]/2.e4, -ge);

  for(int i=0;i<Ee.size();i++)    if (Eeprime[i] < Ee_break && Eeprime[i]>Ee[0])
    {            electrons[i] *= pow( Eeprime[i] / Ee_break, -(ge1-ge));

     if (debug==1)
       cout<<" Ee="<<Ee[i]<< " Phi="<<Phi<<" electron energy + Phi below break="<<Eeprime[i]
           <<" factor="<<pow(Eeprime[i]/Ee_break, -(ge1-ge))
           <<" electrons="<<electrons[i]<<endl;
   }

  // positrons
  double gp = 2.77;
  positrons = 0;
  //  positrons= 2.02e-10* pow(Eeprime /2.e4, -gp);

  for(int i=0;i<Ee.size();i++)    if (Eeprime[i]>Ee[0])   positrons[i]= 2.02e-10* pow(Eeprime[i] /2.e4, -gp);

  for(int i=0;i<Ee.size();i++)    if (Eeprime[i] < Ep_break && Eeprime[i]>Ee[0] )
    {            positrons[i] *= pow( Eeprime[i] / Ep_break, -(gp1-gp));

     if (debug==1)
       cout<<" Ee="<<Ee[i]<< " Phi="<<Phi<<" positron energy + Phi below break="<<Eeprime[i]
           <<" factor="<<pow(Eeprime[i]/Ep_break, -(gp1-gp))
           <<" positrons="<<positrons[i]<<endl;
   }

  result = electrons * electrons_factor + positrons * positrons_factor;


  result *= pow(Ee/(Ee+Phi),2.0);
 }
///////////////////////////////////////////////////////////////

if(name=="imos_electrons_modulated1")
 { 
    
   // electron interstellar spectrum as in   Abdo etal 2011 ApJ 734, 116 (2011) equation (4)
   // modulation as in Abdo etal 2011 ApJ 734, 116 (2011) equations (2) and (3)
   // Phi1: cycles 20/22  Phi2:cycle 21
   // NB   only valid for Phi=400 !
   
   double Phi0    =parameters[0];

   double Phi1,Phi2;
   double rb=100.;
   double r0=10.;

   if(r>=r0)Phi1=Phi0/1.88*          (pow(r,-0.4)-pow(rb,-0.4)); 
   if(r< r0)Phi1=Phi0/1.88*(0.24+8.0*(pow(r,-0.1)-pow(r0,-0.1)));

   double Phi=Phi1; // since using spectrum outside heliosphere

  if(debug==1)cout<<"CR_Spectrum name="<<name<<" r="<<r<<" Phi="<<Phi<<endl;

  double a =160.24;
  double b =  7.0 ;
  double c = -1.2 ;
  double d =  3.03;
  
  valarray<double> EeGeVprime=(Ee+Phi)/1e3;
                                                    result    =  a * pow(b+c, -d) * pow(EeGeVprime/b,    -3.0);
 for(int i=0;i<Ee.size();i++) if(EeGeVprime[i] >=b) result[i] =  a *                pow(EeGeVprime[i]+c, -d  );

 
 

  result*=1.0e-7; // m^-2 GeV^-1 -> cm-2 MeV-1
  result *= pow(Ee/(Ee+Phi),2.0);

 if(debug==1)
   for(int i=0;i<Ee.size();i++) cout<<"imos model 1 r="<<r
   <<" Ee="<<Ee[i]<<" result="<<  result[i]<<" "<<result[i]*1e7*pow(Ee[i]/1e3,3.0)<<endl;

 }
///////////////////////////////////////////////////////////////

if(name=="imos_electrons_modulated2")
 { 
    
   // electron interstellar spectrum as in   Abdo etal 2011 ApJ 734, 116 (2011) equation (4)
   // modulation as in Abdo etal 2011 ApJ 734, 116 (2011) equations (2) and (3)
   // Phi1: cycles 20/22  Phi2:cycle 21
   // NB   only valid for Phi=400 !
   
   double Phi0    =parameters[0];

   double Phi1,Phi2;
   double rb=100.;
   double r0=10.;

   Phi2 = Phi0 * (pow(r,-0.1)-pow(rb,-0.1))/(1.0-pow(rb,-0.1));

   double Phi=Phi2; // since using spectrum outside heliosphere

  if(debug==1)cout<<"CR_Spectrum name="<<name<<" r="<<r<<" Phi="<<Phi<<endl;

  double a =160.24;
  double b =  7.0 ;
  double c = -1.2 ;
  double d =  3.03;
  
 

 valarray<double> EeGeVprime=(Ee+Phi)/1e3;
                                                    result    =  a * pow(b+c, -d) * pow(EeGeVprime/b,    -3.0);
 for(int i=0;i<Ee.size();i++) if(EeGeVprime[i] >=b) result[i] =  a *                pow(EeGeVprime[i]+c, -d  );


 

  result*=1.0e-7; // m^-2 GeV^-1 -> cm-2 MeV-1
  result *= pow(Ee/(Ee+Phi),2.0);

 if(debug==1)
   for(int i=0;i<Ee.size();i++) cout<<"imos model 2 r="<<r
   <<" Ee="<<Ee[i]<<" result="<<  result[i]<<" "<<result[i]*1e7*pow(Ee[i]/1e3,3.0)<<endl;
 }

///////////////////////////////////////////////////////////////

if(name=="imos_electrons_modulated3")
 { 
    
   // electron interstellar spectrum as in   Abdo etal 2011 ApJ 734, 116 (2011) equation (4)
   // modulation as in Abdo etal 2011 ApJ 734, 116 (2011) equations (2) and (3)
   // Phi1: cycles 20/22  Phi2:cycle 21
   // model 3: Phi1 for r> 1 AU, constant r<1 AU
   // NB   only valid for Phi=400 !
 
   
   double Phi0    =parameters[0];

   double Phi1,Phi2;
   double rb=100.;
   double r0=10.;

   
   if(r>=r0)Phi1=Phi0/1.88*          (pow(r,-0.4)-pow(rb,-0.4)); 
   if(r< r0)Phi1=Phi0;

   double Phi=Phi1; // since using spectrum outside heliosphere
   

  if(debug==1)cout<<"CR_Spectrum name="<<name<<" r="<<r<<" Phi="<<Phi<<endl;

  double a =160.24;
  double b =  7.0 ;
  double c = -1.2 ;
  double d =  3.03;

  
  valarray<double> EeGeVprime=(Ee+Phi)/1e3;
                                                    result    =  a * pow(b+c, -d) * pow(EeGeVprime/b,    -3.0);
 for(int i=0;i<Ee.size();i++) if(EeGeVprime[i] >=b) result[i] =  a *                pow(EeGeVprime[i]+c, -d  );

 

  result *= 1.0e-7;                 // m^-2 GeV^-1 -> cm-2 MeV-1
  result *= pow(Ee/(Ee+Phi),2.0);   // modulation

  if(debug==1)
   for(int i=0;i<Ee.size();i++) cout<<"imos model 3 r="<<r
   <<" Ee="<<Ee[i]<<" result="<<  result[i]<<" "<<result[i]*1e7*pow(Ee[i]/1e3,3.0)<<endl;


 }
///////////////////////////////////////////////////////////////////

 if(name=="Fermi2011_linear")
 { // Ackermann etal Fermi-LAT arXiv:1109.052 electrons fit, converting from GeV m^2 to MeV cm^2 
  double ge=3.19;
  result = 2.07e-9* pow(Ee/2.e4, -ge);
  result *= r; // r in AU so linear increase
 }
 /////////////////////////////////////////////////////////////////

  if(debug==1)
    for(int i=0;i<Ee.size();i++) cout<<"CR_Spectrum: name="<<name<<" Phi0="<<parameters[0]<<" r="<<r<<" Ee= "<<Ee[i]<<" result="<<result[i]<<endl;

  if(debug==1)  cout<<"<<CR_Spectrum::fluxe"<<endl;

  return result;
}
//////////////////////////////

void CR_Spectrum::test()
{
  cout<< "CR_Spectrum::test"<<endl;

  CR_Spectrum leptonspectrum;
  double Ee_min=1.0;
  double Ee_max=1e6;
  double Ee_factor=pow(10.,0.5);

  int     n_grid = int(log(Ee_max/Ee_min)/log(Ee_factor)+ 1.001) ;
  valarray<double> Ee;
  valarray<double> fluxe;
  Ee.resize(n_grid);
  fluxe.resize(n_grid);

  for(int i=0;i<Ee.size();i++) Ee[i]= Ee_min*pow(Ee_factor,i);
  cout<<"electron grid size="<<Ee.size()<<endl;
  cout<<"Ee= ";    for(int i=0;i<Ee.size();i++)    cout<<    Ee[i]<<" " ; cout<<endl;

  string name="powerlaw";
  valarray<double> parameters;
  parameters.resize(10);
  
  double constant=1.0;
  double ge      =3.0;
  double Phi0    =500.;

  parameters[0]=Phi0;
  parameters[1]=constant;
  parameters[2]=ge;

  int debug=1;
  leptonspectrum.init(name,parameters,debug);

  double r=0.5;
  fluxe=leptonspectrum.fluxe(Ee,r);

  for(int i=0;i<Ee.size();i++) cout<<  "Ee= "<<  Ee[i]<<" fluxe="<<fluxe[i]<<endl;


  name="Fermi2011_electrons_positrons_modulated1";
  leptonspectrum.init(name,parameters,debug);

   r=0.5;
  fluxe=leptonspectrum.fluxe(Ee,r);

  for(int i=0;i<Ee.size();i++) cout<<  "Ee= "<<  Ee[i]<<" fluxe="<<fluxe[i]<<endl;

  name="Fermi2011_electrons_positrons";
  parameters[0]=Phi0;
  parameters[1]=3.e3;// Ee_break
  parameters[2]=1.0; // ge1=electron index below break
  parameters[3]=5.e3;// Ep_break
  parameters[4]=2.0; // gp1=positron index below break

  parameters[5]=1.0;//  electrons factor
  parameters[6]=1.0;//  positrons factor


  leptonspectrum.init(name,parameters,debug);

   r=1.0;
  fluxe=leptonspectrum.fluxe(Ee,r);

  for(int i=0;i<Ee.size();i++) cout<<  "Ee= "<<  Ee[i]<<" fluxe="<<fluxe[i]<<endl;

  name="Fermi2011_electrons_positrons_modulated1";
  leptonspectrum.init(name,parameters,debug);

  
  fluxe=leptonspectrum.fluxe(Ee,r);

  for(int i=0;i<Ee.size();i++) cout<<  "Ee= "<<  Ee[i]<<" fluxe="<<fluxe[i]<<endl;

  name="Fermi2011_electrons_positrons_modulated2";
  leptonspectrum.init(name,parameters,debug);

  fluxe=leptonspectrum.fluxe(Ee,r);

  for(int i=0;i<Ee.size();i++) cout<<  "Ee= "<<  Ee[i]<<" fluxe="<<fluxe[i]<<endl;

  r=2.0;
  name="Fermi2011_electrons_positrons";
  leptonspectrum.init(name,parameters,debug);
  fluxe=leptonspectrum.fluxe(Ee,r);

  for(int i=0;i<Ee.size();i++) cout<<  "Ee= "<<  Ee[i]<<" fluxe="<<fluxe[i]<<endl;

  name="Fermi2011_electrons_positrons_modulated1";
  leptonspectrum.init(name,parameters,debug);
  fluxe=leptonspectrum.fluxe(Ee,r);

  for(int i=0;i<Ee.size();i++) cout<<  "Ee= "<<  Ee[i]<<" fluxe="<<fluxe[i]<<endl;

  name="Fermi2011_electrons_positrons_modulated2";
  leptonspectrum.init(name,parameters,debug);
  fluxe=leptonspectrum.fluxe(Ee,r);

  for(int i=0;i<Ee.size();i++) cout<<  "Ee= "<<  Ee[i]<<" fluxe="<<fluxe[i]<<endl;




  // CR flux without r

  double E_min=1.0; // GeV
  double E_max=1e6;
  double E_factor=pow(10.,0.2);

  n_grid = int(log(E_max/E_min)/log(E_factor)+ 1.001) ;
  valarray<double> E;
  valarray<double> flux;
  E    .resize(n_grid);
  flux .resize(n_grid);

  for(int i=0;i<E.size();i++) E[i]= E_min*pow(E_factor,i);
  cout<<"energy grid size="<<E.size()<<endl;
  cout<<"E= ";    for(int i=0;i<E.size();i++)    cout<<    E[i]<<" " ; cout<<endl;

  valarray<double> EMeV; // to test MeV versions
  EMeV .resize(n_grid);
  EMeV = E*1.e3; // GeV->MeV

  name="powerlaw";
  
  parameters.resize(10);
  
  constant=1.0;
  ge      =3.0;
  double Phi    =500.; // not Phi0 as above !

  parameters[0]=Phi;
  parameters[1]=constant;
  parameters[2]=ge;

  debug=1;
  CR_Spectrum cr_spectrum;
  cr_spectrum.init(name,parameters,debug);

  flux=cr_spectrum.flux(E);
  for(int i=0;i<E.size();i++) cout<<"CR_Spectrum.test: name="<<name<<  " E= "<<  E[i]<<" flux="<<flux[i]<<endl;


  name="Shikaze protons";
   
  cr_spectrum.init(name,parameters,debug);
  flux=cr_spectrum.flux(E);
  for(int i=0;i<E.size();i++) cout<<"CR_Spectrum.test: name="<<name<<  " E= "<<  E[i]<<" flux="<<flux[i]<<endl;

  name="Shikaze protons MeV";
   
  cr_spectrum.init(name,parameters,debug);
  flux=cr_spectrum.flux(EMeV);
  for(int i=0;i<E.size();i++) cout<<"CR_Spectrum.test: name="<<name<<  " EMeV= "<<  EMeV[i]<<" flux="<<flux[i]<<endl;


  name="Shikaze Helium";
   
  cr_spectrum.init(name,parameters,debug);
  flux=cr_spectrum.flux(E);
  for(int i=0;i<E.size();i++) cout<<"CR_Spectrum.test: name="<<name<<  " E= "<<  E[i]<<" flux="<<flux[i]<<endl;

  name="Shikaze Helium MeV";
   
  cr_spectrum.init(name,parameters,debug);
  flux=cr_spectrum.flux(EMeV);
  for(int i=0;i<E.size();i++) cout<<"CR_Spectrum.test: name="<<name<<  " EMeV= "<<  EMeV[i]<<" flux="<<flux[i]<<endl;



  name="Gaisser protons";
   
  cr_spectrum.init(name,parameters,debug);
  flux=cr_spectrum.flux(E);
  for(int i=0;i<E.size();i++) cout<<"CR_Spectrum.test: name="<<name<<  " E= "<<  E[i]<<" flux="<<flux[i]<<endl;

  name="Gaisser Helium high";
   
  cr_spectrum.init(name,parameters,debug);
  flux=cr_spectrum.flux(E);
  for(int i=0;i<E.size();i++) cout<<"CR_Spectrum.test: name="<<name<<  " E= "<<  E[i]<<" flux="<<flux[i]<<endl;

  name="Gaisser Helium low";
   
  cr_spectrum.init(name,parameters,debug);
  flux=cr_spectrum.flux(E);
  for(int i=0;i<E.size();i++) cout<<"CR_Spectrum.test: name="<<name<<  " E= "<<  E[i]<<" flux="<<flux[i]<<endl;


 return;
}
