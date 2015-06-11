  // Uses gamma-ray (and eventually also other products) hadronic production matrices from
  // Huang, C.-Y.; Park, S.-E.; Pohl, M.; Daniels, C. D. 2007 Astroparticle Physics, Volume 27, Issue 5, p. 429-439.
  // http://cdsads.u-strasbg.fr/cgi-bin/bib_query?2007APh....27..429H
  // More detail in Delahaye et al. 2011 http://adsabs.harvard.edu/abs/2011A%26A...531A..37D

  // Data can be downloaded from http://www.app.physik.uni-potsdam.de/gamma-prod.html
  // data files contain matrices with 201 lines of 374 entries, each entry is a 15 byte number in format .1234567890D-01
  // each of the 201 lines is for one gamma-ray (or other product) energy, 
  // with the entries corresponding to the 374 proton or Helium energies
  // energies are total in GeV/nucleon for proton and Helium, kinetic energy in GeV for products.

  // Implemented here so far:
  // CR Protons + ISM -> gamma rays
  // CR Helium  + ISM -> gamma rays

  // other secondaries will follow later 

// Requires the following files:
// source:
// CR_Spectrum.cc  CR_Spectrum.h  ProductionMatrices.cc  ProductionMatrices.h  test_ProductionMatrices.cc
// data:
// espectra_gamma.decay.he.matrix.data  espectra_gamma.decay.p.matrix.data  prodxsection.he.matrix.data  prodxsection.p.matrix.data


// Uses gamma-ray production function from Kachelriess & Ostapchenko http://arxiv.org/abs/1206.4705
// additional source for  Kachelriess & Ostapchenko 2012:
// frag.f
// additional data   for  Kachelriess & Ostapchenko 2012:
// apfrag.dat  gamfrag.dat
// additional data   for  Kachelriess & Ostapchenko 2013 ppfrag v02 revised April 2013
// gam_ppR  gam_pHe_nm  gam_Hep_nm  gam_HeHeR

// Uses gamma-ray production function from  Dermer 1986 AA 157, 223. Code from C. Dermer modified for c++ by A. Strong
// DERMER_KAMAE_PION_GAMMA_3.cc

 // Sample program is in test_ProductionMatrices.cc
//  to compile and run sample program:
//  Version 1.0: g++ *cc -O3 (also tested with intel icpc)
//  Version 2.0: g++ *cc frag.f -lgfortran -lgfortranbegin  -O3
//  Version 3.0: g++ *cc frag.f *.f90 -lgfortran -lgfortranbegin  -O3
//  Version 4.0:  g++ *cc frag.f *.f90 ../EnergyDispersion/EnergyDispersion.cc -I../EnergyDispersion -I$FITSIO_DIRECTORY/include -L$FITSIO_DIRECTORY/lib -lgfortran -lgfortranbegin -lcfitsio  -O3
//  Version 5.0:  g++ CR_Spectrum.cc DERMER_KAMAE_PION_GAMMA_13.cc kamae.cc ProductionMatrices.cc test_ProductionMatrices.cc frag.f *.f90 ../EnergyDispersion/EnergyDispersion.cc -I../EnergyDispersion -I$FITSIO_DIRECTORY/include -L$FITSIO_DIRECTORY/lib -lgfortran -lgfortranbegin -lcfitsio

//  needs also e.g. setenv FITSIO_DIRECTORY /afs/ipp-garching.mpg.de/home/a/aws/propagate/c/cfitsio/3.26/gcc_sles11_olga2/cfitsio

//  ./a.out  (lots of screen output)
//  to see just the sample output:
// ./a.out | grep sample

// Version 1.0  13 Apr 2012 Huang matrices
// Version 2.0  25 Jul 2012 Kachelriess & Ostapchenko functions added
// Version 3.0  05 Feb 2013 Kachelriess & Ostapchenko functions version 2 added. NB *.f90 copies of f95
// Version 5.0  07 May 2013 new Dermer cross-sections ( DERMER_KAMAE_PION_GAMMA_13.cc replaces  DERMER_KAMAE_PION_GAMMA_3.cc)
// Andy Strong, MPE, Garching, aws@mpe.mpg.de
// Available from http://www.mpe.mpg.de/~aws/propagate.html
 

//////////////////////////////////////////////////////////////////////////////////////////////////

using namespace std;
#include<iostream>
#include<fstream>
#include<cstring>
#include<cstdlib>
#include<vector>
#include<valarray>

#include"EnergyDispersion.h"
#include"ProductionMatrices.h"

#include"CR_Spectrum.h"

// for kamae.cc
valarray<double>     kamae_gamma_param_nd   (double Tp);
valarray<double>     kamae_gamma_param_diff (double Tp);
valarray<double>     kamae_gamma_param_delta(double Tp);
valarray<double>     kamae_gamma_param_res  (double Tp);

double  kamae_nd   (double E,double Tp,valarray<double>param_a,int particle);
double  kamae_diff (double E,double Tp,valarray<double>param_b);
double  kamae_delta(double E,double Tp,valarray<double>param_c);
double  kamae_res  (double E,double Tp,valarray<double>param_d);

// Dermer function in  DERMER_KAMAE_PION_GAMMA_3.cc
double EdsigDermerdE(double Esec, double TP); //Esec, Tp in GeV, Units of mb



/////////////////////////////////////////////////
int  ProductionMatrices::init(string energy_type_,string energy_units_, int debug_, string model_name_,int energy_dispersion_)
{
  cout<<">>ProductionMatrices::init";
  energy_type      = energy_type_;
  energy_units     = energy_units_;
  model_name       = model_name_;
  energy_dispersion=energy_dispersion_;

  cout<<" energy_type      = " <<energy_type;
  cout<<" energy_units     = " <<energy_units<<endl;
  cout<<" model_name       = " <<model_name  <<endl;
  cout<<" energy_dispersion= " <<energy_dispersion  <<endl;

  if (energy_type !="total_energy"&&energy_type !="kinetic_energy"){ cout<<"invalid energy type !"<<endl; exit(1);}
  if (energy_units!="GeV"         &&energy_units!="MeV"           ){ cout<<"invalid energy units!"<<endl; exit(1);}
  if (  model_name  !="Huang"       
      &&model_name  !="Kachelriess" 
      &&model_name  !="Kachelriess2"  
      &&model_name  !="Kachelriess2013"
      &&model_name  !="Kamae_Kachelriess"
      &&model_name  !="Dermer_Kachelriess" 
      &&model_name  !="Dermer_Kachelriess2013" 
      &&model_name  !="Dermer_Kachelriess2013_heavy" 
      &&model_name  !="Kamae"               
      &&model_name  !="Dermer"     
      &&model_name  !="FLUKA"           
     ) { cout<<"invalid model name  !"<<endl; exit(1);}
  
  debug = debug_;

  cout<<"<<ProductionMatrices::init"<<endl;
  return 0;
}
/////////////////////////////////////////////////////

int ProductionMatrices::read(string input_directory)
{



  cout<<" ProductionMatrices::read"<<endl;
  cout<<" debug="<<debug<<endl;

   
  string infile_string;

  int process;

  // read the gamma-ray spectra matrices

  if (model_name=="Huang") //AWS20131104
  {  
   cout<< "reading Huang matrices"<<endl;

 for(process=1;process<=2;process++)
 {

  

  infile_string=input_directory + "/";  
  if(process==1)infile_string+="espectra_gamma.decay.p.matrix.data";
  if(process==2)infile_string+="espectra_gamma.decay.he.matrix.data";

  cout<<"reading from "<<infile_string<<endl;

  ifstream in;
  in.open(infile_string.c_str());

  char *buffer;



  nE_proton=374;
  nE_gamma =201;

  int number_size=15;

  double *values;
  values = new double[nE_proton];

  vector<valarray<double> > matrix;
  matrix.resize(nE_proton);
  for(int i=0;i<nE_proton;i++) matrix[i].resize(nE_gamma);

 // see Schildt c++ Complete Reference page 554: note the +1 needed for getline

  int linelength=nE_proton*number_size;
  buffer =new char[linelength];


 for(int iE_gamma=0;iE_gamma<nE_gamma;iE_gamma++)
 {

  in.getline(buffer,linelength+1);
 
  if(debug==1)cout<<"input buffer="<<buffer<<endl;

  // change number format from fortran D to c++ e

  for(int i=0;i<nE_proton;i++)
  {
  //cout<<"strchr:"<<strchr(buffer,'D')<<endl; //"D" does not compile, 'D' does

   strncpy(strchr(buffer,'D'),"e",1); // copies one character

  }

  if(debug==1)  cout<<"modified buffer="<<buffer<<endl;

  for(int i=0;i<nE_proton;i++)
  {  
   sscanf(&buffer[i*number_size],"%le",&values[i]); // le = double, Schildt p.209
  }

  if(debug==1)
  {
   cout<<"values=";
   for(int i=0;i<nE_proton;i++)cout<<values[i]<<" "; cout<<endl;
  }


  for(int iE_proton=0;iE_proton<nE_proton;iE_proton++)
   matrix[iE_proton][iE_gamma] =   values[iE_proton];

    }//iE_gamma

 in.close();


 E_proton.resize(nE_proton);
 E_gamma .resize(nE_gamma );

 // equation 1 of Huang manual note-1.pdf
 // TOTAL energy in GeV/nucleon
  for(int iE_proton=0;iE_proton<nE_proton;iE_proton++)
    E_proton[iE_proton] = 1.24 * pow( (1.0+0.05), iE_proton-1);

 // equation 2 of Huang manual note-1.pdf
  for(int iE_gamma =0;iE_gamma <nE_gamma ;iE_gamma ++)
    E_gamma[iE_gamma]= exp(log(0.01) + iE_gamma*(log(1.e8)-log(0.01))/nE_gamma);


 if(debug==2)
 {

  cout<<"proton total energies (GeV/nucleon)="<<endl;
  for(int iE_proton=0;iE_proton<nE_proton;iE_proton++)cout<< E_proton[iE_proton]<<" "; cout<<endl;

  cout<<"gamma  energies (GeV)="<<endl;
  for(int iE_gamma =0;iE_gamma <nE_gamma ;iE_gamma ++)cout<< E_gamma [iE_gamma ]<<" "; cout<<endl;




  cout<<"matrix="<<endl;

  for(int iE_gamma =0;iE_gamma <nE_gamma ;iE_gamma ++)
  {
   cout<<"iE_gamma="<<iE_gamma<<endl;

  for(int iE_proton=0;iE_proton<nE_proton;iE_proton++)
         cout<<matrix[iE_proton][iE_gamma]<<" ";
    
  cout<<endl;
  }

 }//debug==2


 if(process==1) gamma_decay_p_matrix =matrix;
 if(process==2) gamma_decay_he_matrix=matrix;


 }// process


 ///////////////////////////////////////////////////////////////////////

  // read the total inelastic cross-sections

   

 prodxsection_p_matrix .resize(nE_proton);
 prodxsection_he_matrix.resize(nE_proton);

 for(process=1;process<=2;process++)
 {
  infile_string=input_directory + "/";  
  if(process==1)infile_string+="prodxsection.p.matrix.data";
  if(process==2)infile_string+="prodxsection.he.matrix.data";

 

  cout<<"reading from "<<infile_string<<endl;

  ifstream in;
  in.open(infile_string.c_str());

  char *buffer;
  

  int number_size=15;

  double *values;
  values = new double[nE_proton];

  vector<valarray<double> > matrix;
  matrix.resize(nE_proton);
  for(int i=0;i<nE_proton;i++) matrix[i].resize(1);

 // see Schildt c++ Complete Reference page 554: note the +1 needed for getline

  int linelength=nE_proton*number_size;
  buffer =new char[linelength];

  in.getline(buffer,linelength+1);
 
  if(debug==1)cout<<"prodxsection input buffer="<<buffer<<endl;

  // change number format from fortran D to c++ e

  for(int i=0;i<nE_proton;i++)
  {
  //cout<<"strchr:"<<strchr(buffer,'D')<<endl; //"D" does not compile, 'D' does

   strncpy(strchr(buffer,'D'),"e",1); // copies one character

  }

  if(debug==1)  cout<<" prodx modified buffer="<<buffer<<endl;

  for(int i=0;i<nE_proton;i++)
  {  
   sscanf(&buffer[i*number_size],"%le",&values[i]); // le = double, Schildt p.209
  }

  if(debug==1)
  {
   cout<<"prodx values=";
   for(int i=0;i<nE_proton;i++)cout<<values[i]<<" "; cout<<endl;
  }


 in.close();


 if(debug==2)
 {

  cout<<"proton total energies (GeV/nucleon)="<<endl;
  for(int iE_proton=0;iE_proton<nE_proton;iE_proton++)cout<< E_proton[iE_proton]<<" "; cout<<endl;

  cout<<"prodx values="<<endl;

  
  for(int iE_proton=0;iE_proton<nE_proton;iE_proton++)  cout<<values[iE_proton]<<" ";
    
  cout<<endl;
  

 }//debug==2


 if(process==1) for(int iE_proton=0;iE_proton<nE_proton;iE_proton++) prodxsection_p_matrix [iE_proton]=values[iE_proton];
 if(process==2) for(int iE_proton=0;iE_proton<nE_proton;iE_proton++) prodxsection_he_matrix[iE_proton]=values[iE_proton];


 }// process



  } // model_name=="Huang"



 // FLUKA matrices

if(model_name=="FLUKA")//AWS20131104
{
 cout<< "reading FLUKA matrices"<<endl;

 for(process=1;process<=2;process++)
 {

  
  infile_string=input_directory + "/";  
  if(process==1)infile_string+="fluka_matrix_aws_1.dat";
  if(process==2)infile_string+="fluka_matrix_aws_1.dat"; // just for consistency for now

  cout<<"reading from "<<infile_string<<endl;

  ifstream in;
  in.open(infile_string.c_str());
  in >> nE_proton_FLUKA >> nE_gamma_FLUKA;
  cout<<" nE_proton_FLUKA="<< nE_proton_FLUKA<<"  nE_gamma_FLUKA="<< nE_gamma_FLUKA<<endl;

  E_proton_FLUKA       .resize( nE_proton_FLUKA);
  E_gamma_FLUKA        .resize( nE_gamma_FLUKA );

  gamma_p_FLUKA_matrix .resize(nE_proton_FLUKA);
  gamma_he_FLUKA_matrix.resize(nE_proton_FLUKA); // just for consistency for now

  for(int iE_proton=0;iE_proton<nE_proton_FLUKA;iE_proton++)  gamma_p_FLUKA_matrix [iE_proton].resize(nE_gamma_FLUKA);
  for(int iE_proton=0;iE_proton<nE_proton_FLUKA;iE_proton++)  gamma_he_FLUKA_matrix[iE_proton].resize(nE_gamma_FLUKA); // just for consistency for now


  for(int iE_proton=0;iE_proton<nE_proton_FLUKA;iE_proton++) in>>E_proton_FLUKA[iE_proton];
  for(int iE_gamma =0;iE_gamma <nE_gamma_FLUKA ;iE_gamma++ ) in>>E_gamma_FLUKA [iE_gamma];

  cout<<"E_proton_FLUKA = ";
  for(int iE_proton=0;iE_proton<nE_proton_FLUKA;iE_proton++)  cout<<E_proton_FLUKA[iE_proton]<<" "; cout<<endl;
  cout<<"E_gamma_FLUKA = ";
  for(int iE_gamma=0;iE_gamma<nE_gamma_FLUKA ;iE_gamma++ ) cout<<E_gamma_FLUKA [iE_gamma]<<" "; cout<<endl;

  for(int iE_proton=0;iE_proton<nE_proton_FLUKA;iE_proton++)
  for(int iE_gamma =0;iE_gamma <nE_gamma_FLUKA ;iE_gamma++ ) in>>gamma_p_FLUKA_matrix[iE_proton][iE_gamma];



  gamma_he_FLUKA_matrix = gamma_p_FLUKA_matrix; // just for consistency for now

 
  cout<<"FLUKA matrix="<<endl;

  for(int iE_gamma =0;iE_gamma <nE_gamma_FLUKA ;iE_gamma ++)
  {
    cout<<"FLUKA matrix iE_gamma="<<iE_gamma<<" E_gamma_FLUKA="<<E_gamma_FLUKA [iE_gamma]<<": " ;

  for(int iE_proton=0;iE_proton<nE_proton_FLUKA;iE_proton++)
         cout<<gamma_p_FLUKA_matrix[iE_proton][iE_gamma]<<" ";
    
  cout<<endl;
  }

 }

 } // model_name=="FLUKA"

 matrices_read=1;

  return 0;
};

/////////////////////////////////////////////////////

int ProductionMatrices::print()
{

  cout<<"ProductionMatrices::print"<<endl;

  // print input matrices

  cout<<"proton total energies (GeV/nucleon)="<<endl;
  for(int iE_proton=0;iE_proton<nE_proton;iE_proton++)cout<< E_proton[iE_proton]<<" "; cout<<endl;

  cout<<"gamma  energies (GeV)="<<endl;
  for(int iE_gamma =0;iE_gamma <nE_gamma ;iE_gamma ++)cout<< E_gamma [iE_gamma ]<<" "; cout<<endl;

  cout<<" input prodxsection_p_matrix  = ";
  for(int iE_proton=0;iE_proton<nE_proton;iE_proton++) cout<< prodxsection_p_matrix  [iE_proton]<<" "; cout<<endl;

  cout<<" input prodxsection_he_matrix  = ";
  for(int iE_proton=0;iE_proton<nE_proton;iE_proton++) cout<< prodxsection_he_matrix [iE_proton]<<" "; cout<<endl;

  cout<<" input gamma_decay_p_matrix="<<endl;

  for(int iE_gamma =0;iE_gamma <nE_gamma ;iE_gamma ++)
  {
   cout<<"iE_gamma="<<iE_gamma<<" E_gamma="<<E_gamma[iE_gamma] <<endl;

  for(int iE_proton=0;iE_proton<nE_proton;iE_proton++)
         cout<<  gamma_decay_p_matrix[iE_proton][iE_gamma]<<" ";
    
  cout<<endl;
  }

  cout<<" input gamma_decay_he_matrix="<<endl;

  for(int iE_gamma =0;iE_gamma <nE_gamma ;iE_gamma ++)
  {
    cout<<"iE_gamma="<<iE_gamma<<" E_gamma="<<E_gamma[iE_gamma]<<endl;

  for(int iE_proton=0;iE_proton<nE_proton;iE_proton++)
         cout<<  gamma_decay_he_matrix[iE_proton][iE_gamma]<<" ";
    
  cout<<endl;
  }

  // print rebinned matrices
  cout<<endl;
  cout<<"rebinned proton total energies ("<<energy_units<<"/nucleon)="<<endl;
  for(int iE_proton_rebin=0;iE_proton_rebin<nE_proton_rebin;iE_proton_rebin++)cout<< E_proton_rebin[iE_proton_rebin]<<" "; cout<<endl;

  cout<<endl;
  cout<<"rebinned gamma  energies ("<<energy_units<<")="<<endl;
  for(int iE_gamma_rebin =0;iE_gamma_rebin <nE_gamma_rebin ;iE_gamma_rebin ++)cout<< E_gamma_rebin [iE_gamma_rebin ]<<" "; cout<<endl;

  cout<<endl;
  cout<<" view gamma_decay_p_matrix_rebin="<<endl;

  for(int iE_gamma_rebin =0;iE_gamma_rebin <nE_gamma_rebin ;iE_gamma_rebin ++)
  {
    cout<<"iE_gamma_rebin="<<iE_gamma_rebin<<" E_gamma_rebin="<<E_gamma_rebin[iE_gamma_rebin]<<" "<<energy_units<<endl;

  for(int iE_proton_rebin=0;iE_proton_rebin<nE_proton_rebin;iE_proton_rebin++)
         cout<<  gamma_decay_p_matrix_rebin[iE_proton_rebin][iE_gamma_rebin]<<" ";
    
  cout<<endl;
  }
  cout<<endl;
  cout<<" view gamma_decay_he_matrix_rebin="<<endl;

  for(int iE_gamma_rebin =0;iE_gamma_rebin <nE_gamma_rebin ;iE_gamma_rebin ++)
  {
    cout<<"iE_gamma_rebin="<<iE_gamma_rebin<<" E_gamma_rebin="<<E_gamma_rebin[iE_gamma_rebin]<<" "<<energy_units<<endl;

  for(int iE_proton_rebin=0;iE_proton_rebin<nE_proton_rebin;iE_proton_rebin++)
         cout<<  gamma_decay_he_matrix_rebin[iE_proton_rebin][iE_gamma_rebin]<<" ";
    
  cout<<endl;
  }


  // print emissivities
  if(gamma_decay_p_emissivity_rebin.size() >0) // in case not computed
  {
  cout<<"  gamma_decay_p_emissivity_rebin="<<endl;
  cout<<"  energy units="<<energy_units<<endl;

  for(int iE_gamma_rebin =0;iE_gamma_rebin <nE_gamma_rebin ;iE_gamma_rebin ++)
    cout<<"iE_gamma_rebin="<<iE_gamma_rebin<<" E_gamma_rebin="<<E_gamma_rebin[iE_gamma_rebin]
        <<" "<<energy_units<<" proton emissivity= "<< gamma_decay_p_emissivity_rebin[iE_gamma_rebin]
	<<" sr-1 s-1 "<<energy_units<<"-1 H atom-1 "
        << " *E^2: "<< gamma_decay_p_emissivity_rebin[iE_gamma_rebin]*pow(E_gamma_rebin[iE_gamma_rebin],2)*1e3
        <<" MeV sr-1 s-1  H atom-1 "
        <<endl;

  cout<<"  gamma_decay_he_emissivity_rebin="<<endl;

  for(int iE_gamma_rebin =0;iE_gamma_rebin <nE_gamma_rebin ;iE_gamma_rebin ++)
    cout<<"iE_gamma_rebin="<<iE_gamma_rebin<<" E_gamma_rebin="<<E_gamma_rebin[iE_gamma_rebin]
        <<" "<<energy_units<<"  helium emissivity= "<< gamma_decay_he_emissivity_rebin[iE_gamma_rebin]
	<<" sr-1 s-1 "<<energy_units<<"-1 H atom-1 "
        << " *E^2: "<< gamma_decay_he_emissivity_rebin[iE_gamma_rebin]*pow(E_gamma_rebin[iE_gamma_rebin],2)*1e3
        <<" MeV sr-1 s-1  H atom-1 "
        <<endl;

  cout<<"   gamma_decay_p_emissivity_rebin + gamma_decay_he_emissivity_rebin="<<endl;

  for(int iE_gamma_rebin =0;iE_gamma_rebin <nE_gamma_rebin ;iE_gamma_rebin ++)
    cout<<"iE_gamma_rebin="<<iE_gamma_rebin<<" E_gamma_rebin="<<E_gamma_rebin[iE_gamma_rebin]
        <<" GeV proton + helium emissivity= "
        <<  gamma_decay_p_emissivity_rebin[iE_gamma_rebin]+gamma_decay_he_emissivity_rebin[iE_gamma_rebin]
	<<" sr-1 s-1 GeV-1 H atom-1"
        << " *E^2: "
        <<  (gamma_decay_p_emissivity_rebin [iE_gamma_rebin]
	    +gamma_decay_he_emissivity_rebin[iE_gamma_rebin])*pow(E_gamma_rebin[iE_gamma_rebin],2)*1e3
        <<" MeV sr-1 s-1  H atom-1 "
        <<endl;

  }

  return 0;
};

//////////////////////////////////////////////////////////////////

int ProductionMatrices::rebin(valarray<double> E_proton_rebin_, valarray<double> E_gamma_rebin_ )
{
 cout<<">>ProductionMatrices::rebin"<<endl;


 // int debug=0; now in class set by init
 // if(energy_units_=="MeV")debug=1;

 if(matrices_read!= 1) {cout<<" ProductionMatrices:: rebin: matrices have not been not read !"<<endl; exit(1);}

  nE_proton_rebin=E_proton_rebin_.size();
  nE_gamma_rebin = E_gamma_rebin_.size();

  E_proton_rebin.resize(nE_proton_rebin);
  E_gamma_rebin .resize(nE_gamma_rebin );
 
  E_proton_rebin= E_proton_rebin_;
  E_gamma_rebin = E_gamma_rebin_ ;

  ////  energy_units = energy_units_; // now done in init ! AWS20110411

  cout<<"ProductionMatrices::  rebin: energy_units ="<<energy_units<<endl;

  if(energy_units != "GeV" && energy_units != "MeV")
    {cout<<"ProductionMatrices::rebin: energy_units=" <<energy_units<< " must be GeV or MeV !"<<endl; exit(-1);}

  cout<<"ProductionMatrices::rebin nE_proton_rebin="<< nE_proton_rebin<<endl;
  cout<<"ProductionMatrices::rebin nE_gamma_rebin ="<< nE_gamma_rebin <<endl;

  cout<<"ProductionMatrices::rebin E_proton_rebin ["<<energy_units<<"] =";
  for(int iE_proton_rebin=0; iE_proton_rebin<nE_proton_rebin; iE_proton_rebin++)cout<<E_proton_rebin[ iE_proton_rebin]<<" ";
  cout<<endl;

  cout<<"ProductionMatrices::rebin E_gamma_rebin ["<<energy_units<<"] =";
  for(int iE_gamma_rebin=0;  iE_gamma_rebin <nE_gamma_rebin;  iE_gamma_rebin++) cout<<E_gamma_rebin  [ iE_gamma_rebin]<<" ";
  cout<<endl;

  // working input matrix
  vector<valarray<double> > matrix;

  // working input inelastic cross-section
  valarray<double>  prodx;
  if(model_name=="Huang")prodx.resize(nE_proton); //AWS20131104

  //working rebinned matrix
  vector<valarray<double> > matrix_rebin;
  matrix_rebin.resize(nE_proton_rebin);
  for(int i=0;i<nE_proton_rebin;i++) matrix_rebin[i].resize(nE_gamma_rebin);

  gamma_decay_p_matrix_rebin .clear(); // in case called multiple times, for "=" to work
  gamma_decay_he_matrix_rebin.clear(); // in case called multiple times, for "=" to work
  
  double E_gamma_rebin_GeV; // since original matrices are in GeV
  double v;

int process;

  for(process=1;process<=2;process++)
  {
   if(model_name=="Huang")//AWS20131104
   {
    if(process==1) prodx  =  prodxsection_p_matrix;
    if(process==2) prodx  =  prodxsection_he_matrix;

    if(process==1) matrix =  gamma_decay_p_matrix;
    if(process==2) matrix =  gamma_decay_he_matrix;
    }

  for(int i=0;i<nE_proton_rebin;i++) matrix_rebin[i]=0.;
  

  // using formula from Huang etal. note.pdf (see above) in GeV
  //   E_proton[iE_proton] = 1.24 * pow( (1.0+0.05), iE_proton-1);
  // log(E_proton[iE_proton]) = log(1.24)  +log(1.0+0.05)* (iE_proton-1);

  //     E_gamma[iE_gamma]= exp(log(0.01) + iE_gamma*(log(1.e8)-log(0.01))/nE_gamma);
  //log( E_gamma[iE_gamma])=    log(0.01) + iE_gamma*(log(1.e8)-log(0.01))/nE_gamma);

  // the rebinning is either on GeV or MeV scale, while original is GeV, so compute indices accordingly




  for(int iE_proton_rebin=0; iE_proton_rebin<nE_proton_rebin; iE_proton_rebin++)
  {

    double E_proton_total_rebin_GeV; // since original matrices are in GeV and total energy/nucleon

    if(energy_units=="GeV") E_proton_total_rebin_GeV = E_proton_rebin[ iE_proton_rebin];
    if(energy_units=="MeV") E_proton_total_rebin_GeV = E_proton_rebin[ iE_proton_rebin]* 1.e-3; // MeV -> GeV

    double m_p   = 0.938 ;      // proton mass in GeV  (deprecated hard coding of physics constant !)
    if(energy_type=="kinetic_energy")  E_proton_total_rebin_GeV += m_p; // kinetic energy/nucleon-> total energy/nucleon

 
    if ( E_proton_total_rebin_GeV > m_p) // catch case where total energy is less than rest mass
    {

     int iE_proton;


     iE_proton=(log(E_proton_total_rebin_GeV)-log(1.24))/log(1.05) +1; // rebin in GeV or MeV

     if(debug==1) cout<<"E_proton_rebin="<<E_proton_rebin[ iE_proton_rebin]<<" "<<energy_units<<" iE_proton="
		      << iE_proton<<" E_proton="<<E_proton[iE_proton]<<" - "<<E_proton[iE_proton+1] << " GeV"<<endl;

     //    if(iE_proton>=0 && iE_proton<nE_proton-1) moved
     //    { moved
     for(int iE_gamma_rebin =0; iE_gamma_rebin <nE_gamma_rebin ; iE_gamma_rebin++ )
     {

       //      double E_gamma_rebin_GeV; // since original matrices are in GeV moved

      if(energy_units=="GeV") E_gamma_rebin_GeV=E_gamma_rebin[iE_gamma_rebin];
      if(energy_units=="MeV") E_gamma_rebin_GeV=E_gamma_rebin[iE_gamma_rebin]*1.e-3;

      int iE_gamma;
 

       iE_gamma=(log(  E_gamma_rebin_GeV     )- log(0.01))/((log(1.e8)-log(0.01))/nE_gamma); // rebin in GeV or MeV

      if(debug==1) cout<<"E_gamma_rebin="<<E_gamma_rebin[ iE_gamma_rebin]<<" "<<energy_units<<" iE_gamma ="
		       << iE_gamma <<" E_gamma ="<<E_gamma [iE_gamma ]<<" - "<<E_gamma [iE_gamma+1]<<" GeV" <<endl;


      //     double v; moved 

  if(model_name=="Huang")//AWS20131104
  {
    if(iE_proton>=0 && iE_proton<nE_proton-1)
    {
     if(iE_gamma>=0 && iE_gamma<nE_gamma-1)
     {
	double v1=prodx[iE_proton  ] * matrix[iE_proton  ][iE_gamma];
	double v2=prodx[iE_proton+1] * matrix[iE_proton+1][iE_gamma];

	double v3=prodx[iE_proton  ] * matrix[iE_proton  ][iE_gamma+1];
	double v4=prodx[iE_proton+1] * matrix[iE_proton+1][iE_gamma+1];

 
	// first interpolate to E_proton_rebin for each E_gamma
        // differential cross-section (mb GeV-1)  = total inelastic cross-section (mb) * differential spectrum (GeV-1)

        double v12= v1+
                   (v2-v1)
                  *(E_proton_total_rebin_GeV              - E_proton[ iE_proton])           //AWS20120411
                  /(E_proton      [iE_proton+1]           - E_proton[ iE_proton]);



        double v34= v3+
                   (v4-v3)
                  *(E_proton_total_rebin_GeV              - E_proton[ iE_proton])           //AWS20120411
                  /(E_proton      [iE_proton+1]           - E_proton[ iE_proton]);

        
        if(debug==1) cout<<"v1="<<v1<<" v2="<<v2<<" v12="<<v12<<endl;
        if(debug==1) cout<<"v3="<<v3<<" v4="<<v4<<" v34="<<v34<<endl;


	// then interpolate to E_gamma_rebin

               v=  v12+
	          (v34-v12)
                  *(E_gamma_rebin_GeV                   - E_gamma [ iE_gamma ])
                  /(E_gamma      [iE_gamma +1]          - E_gamma [ iE_gamma ]);

   }//if  iE_gamma
  }//if  iE_proton

  } // model_name=="Huang"

     ///////////////////////////////////////////////////////////////////////////////////////


       if(model_name=="Kachelriess")
       {
	 // using original combination of Kamae non-diffractive with QCSJET: underestimated at low energies ?
        double ep,es;
        int id=0;  // gammas
        int reac;  // reaction
        double ethr=50.;// threshold in GeV for Kamae/QCSJET-II. Kamae non-diffractive only (i.e. s_nd+s_delta+s_res)
  
        ep = E_proton_total_rebin_GeV  ;
        es = E_gamma_rebin_GeV;
	double v_ppfrag,v_ppfrag_pHe, v_ppfrag_HeHe;;
        double A1,A2,A_factor;

        double  H_fraction=0.9;// H_fraction defined for Huang, in function emissivity, use the same here
        double He_fraction=1.0-H_fraction;
       

	if(process==1) // CR protons
       {
        reac=0; // pp with Kamae below 50 GeV, QCSJET-II above 50 GeV

        v_ppfrag   = ppfrag_spec_int(ep,es,id,reac);
             

	// approximate correction for He target below 50 GeV
        
        if(ep<ethr) 
	{
	 A1=1.;
         A2=4.;
	 A_factor=pow(A1*A2, 0.8); // Norbury, J. W., & Townsend, L. W. 2007, Nucl. Instrum. Meth. B 254, 187 
         v_ppfrag_pHe = v_ppfrag * A_factor * He_fraction;
	}

        if (ep>=ethr)
	{
         reac=2; // p-He QCSJET-II above 50 GeV

         v_ppfrag_pHe   = ppfrag_spec_int(ep,es,id,reac);
         v_ppfrag_pHe  *= He_fraction;

        }

        v_ppfrag += v_ppfrag_pHe;
	
       }// process==1

  
  

	if(process==2) // CR He
      {


        if(ep<ethr)
	{
         reac=0; // Kamae pp
        
         v_ppfrag     = ppfrag_spec_int(ep,es,id,reac);
         v_ppfrag_HeHe= v_ppfrag;

	 A1=4.;
         A2=1.;
	 A_factor=pow(A1*A2, 0.8); // Norbury, J. W., & Townsend, L. W. 2007, Nucl. Instrum. Meth. B 254, 187 
         v_ppfrag  *= A_factor; // approximate correction for He projectile

	 // approximation for HeHe
	 A1=4.;
         A2=4.;
	 A_factor=pow(A1*A2, 0.8);  
         v_ppfrag_HeHe *= A_factor *He_fraction;

         v_ppfrag += v_ppfrag_HeHe;
	}

        if(ep>=ethr)
	{
         reac=3; // He-p, >=50 GeV only
        
         v_ppfrag   = ppfrag_spec_int(ep,es,id,reac);

         v_ppfrag_HeHe=v_ppfrag;

	 // approximate correction for HeHe, relative to He-p
	 A1=4.;
         A2=4.;
	 A_factor=pow(A1*A2, 0.8)/pow(A1*1.0, 0.8);  
         v_ppfrag_HeHe *= A_factor * He_fraction;

         v_ppfrag += v_ppfrag_HeHe;

	}      

      } // process==2

     // returns dsigma/dlogE = E dsigma/dE, so divide by E. Also for other values for printed comparison

	v_ppfrag       /=es;
	v_ppfrag_pHe   /= es;
	v_ppfrag_HeHe  /= es;

        v_ppfrag     *=             H_fraction; // to correspond to Huang per atom definition assumed in function emissivity
        v_ppfrag_pHe *=             H_fraction;
        v_ppfrag_HeHe*=             H_fraction;


       if(debug==1)
       {
        cout<<"process="<<process;
        if(process==1) cout<<" CR protons ";
        if(process==2) cout<<" CR He      ";
	cout<<" E_proton_total_rebin_GeV=" <<E_proton_total_rebin_GeV
            <<       " E_gamma_rebin_GeV=" <<E_gamma_rebin_GeV
            <<" comparing Huang="<<v<<" with Kachelriess="<< v_ppfrag
            <<" K/H="<<v_ppfrag/v;
	if(process==1) cout<<" v_ppfrag_pHe=" <<v_ppfrag_pHe;
	if(process==2) cout<<" v_ppfrag_HeHe="<<v_ppfrag_HeHe;
        cout<<endl;

        if(process==1)cout<<" CR protons ";
        if(process==2)cout<<" CR He      ";
	cout<<" E_projectile =" <<E_proton_total_rebin_GeV
            <<" E_gamma=" <<E_gamma_rebin_GeV
            <<" comparing Huang="<<v<<" with Kachelriess="<< v_ppfrag
            <<" K/H="<<v_ppfrag/v
            <<endl;
       }//if debug==1


        v= v_ppfrag;
	

     }// model_name=="Kachelriess"


       ////////////////////////////////////////////////////////////
       // Kamae_Kachelriess is the same as "Kachelriess2


       if(model_name=="Kachelriess2" || model_name=="Kamae_Kachelriess" || model_name=="Kamae") 
       {
	 // using combination of   Kamae  with QCSJET
        double ep,es;
        int id=0;  // gammas
        int reac;  // reaction
        double ethr=20.;// threshold in GeV used here for Kamae/QCSJET-II. Should be parameter. 
        if( model_name=="Kamae") ethr=1e6;// use Kamae only so set threshold arbitrarily high 
  
        ep = E_proton_total_rebin_GeV  ;
        es = E_gamma_rebin_GeV;

	double v_ppfrag,v_ppfrag_pHe, v_ppfrag_HeHe;;
        double A1,A2,A_factor;

        double  H_fraction=0.9;// H_fraction defined for Huang, in function emissivity, use the same here
        double He_fraction=1.0-H_fraction;
       
        valarray<double> param_a; param_a.resize(9);
        valarray<double> param_b; param_b.resize(8);
        valarray<double> param_c; param_c.resize(5);
        valarray<double> param_d; param_d.resize(5);


       if(process==1) // CR protons
       {

           
        
        if(ep<ethr) 
	{
	  // Kamae function uses kinetic energy in GeV
          // returns cross-section in mb/GeV times gamma-ray energy in GeV
	
         double Tp=ep-m_p; // Kamae notation for kinetic energy of proton in GeV
         param_a=kamae_gamma_param_nd   (Tp);
         param_b=kamae_gamma_param_diff (Tp);
         param_c=kamae_gamma_param_delta(Tp);
         param_d=kamae_gamma_param_res  (Tp);

         int particle=0; // gammas
         double s_nd    = kamae_nd   (es,Tp,param_a,particle);
         double s_diff  = kamae_diff (es,Tp,param_b);
         double s_delta = kamae_delta(es,Tp,param_c);
         double s_res   = kamae_res  (es,Tp,param_d);
   

        // use this notation for consistency with rest of Kachelriess section. 
         v_ppfrag= (s_nd + s_diff + s_delta + s_res);

         if(debug==1)
         cout<<"Kamae routines: CR p  Tp="<<Tp<<" es="<<es<<" s_nd="<<s_nd<<" s_diff="<<s_diff<<" s_delta="<<s_delta<<" s_res="<<s_res<<" v_ppfrag="<<v_ppfrag<<endl;


	 A1=1.;
         A2=4.;
	 A_factor=pow(A1*A2, 0.8); // Norbury, J. W., & Townsend, L. W. 2007, Nucl. Instrum. Meth. B 254, 187 
         v_ppfrag_pHe = v_ppfrag * A_factor * He_fraction;

	}

        if (ep>=ethr)
	{
         reac=1; // pp QCSJET-II, >10 GeV only

         v_ppfrag   = ppfrag_spec_int(ep,es,id,reac);

         reac=2; // p-He QCSJET-II >10 GeV only

         v_ppfrag_pHe   = ppfrag_spec_int(ep,es,id,reac);
         v_ppfrag_pHe  *= He_fraction;


        }

          v_ppfrag += v_ppfrag_pHe;
	
       }// process==1

  
  

	if(process==2) // CR He
      {


        if(ep<ethr)
	{
         double Tp=ep-m_p; // Kamae notation for kinetic energy of proton in GeV

         param_a=kamae_gamma_param_nd   (Tp);
         param_b=kamae_gamma_param_diff (Tp);
         param_c=kamae_gamma_param_delta(Tp);
         param_d=kamae_gamma_param_res  (Tp);

         int particle=0; // gammas
         double s_nd    = kamae_nd   (es,Tp,param_a,particle);
         double s_diff  = kamae_diff (es,Tp,param_b);
         double s_delta = kamae_delta(es,Tp,param_c);
         double s_res   = kamae_res  (es,Tp,param_d);
   

        // use this notation for consistency with rest of Kachelriess section
         v_ppfrag= (s_nd + s_diff + s_delta + s_res);

         if(debug==1)
         cout<<"Kamae routines: CR He: Tp="<<Tp<<" es="<<es<<" s_nd="<<s_nd<<" s_diff="<<s_diff<<" s_delta="<<s_delta<<" s_res="<<s_res<<" v_ppfrag="<<v_ppfrag<<endl;

         v_ppfrag_HeHe= v_ppfrag;

	 A1=4.;
         A2=1.;
	 A_factor=pow(A1*A2, 0.8); // Norbury, J. W., & Townsend, L. W. 2007, Nucl. Instrum. Meth. B 254, 187 
         v_ppfrag  *= A_factor; // approximate correction for He projectile

	 // approximation for HeHe
	 A1=4.;
         A2=4.;
	 A_factor=pow(A1*A2, 0.8);  
         v_ppfrag_HeHe *= A_factor *He_fraction;

         v_ppfrag += v_ppfrag_HeHe;
	}

        if(ep>=ethr)
	{
         reac=3; // He-p, >=10 GeV only
        
         v_ppfrag   = ppfrag_spec_int(ep,es,id,reac);

         v_ppfrag_HeHe=v_ppfrag;

	 // approximate correction for HeHe, relative to He-p
	 A1=4.;
         A2=4.;
	 A_factor=pow(A1*A2, 0.8)/pow(A1*1.0, 0.8);  
         v_ppfrag_HeHe *= A_factor * He_fraction;

         v_ppfrag += v_ppfrag_HeHe;

	}      

      } // process==2

     // returns dsigma/dlogE = E dsigma/dE, so divide by E. Also for other values for printed comparison

	v_ppfrag       /= es;
	v_ppfrag_pHe   /= es;
	v_ppfrag_HeHe  /= es;

        v_ppfrag     *=             H_fraction; // to correspond to Huang per atom definition assumed in function emissivity: but so it is not per atom ! 
        v_ppfrag_pHe *=             H_fraction; // but cancels in emissivity where division by H_fraction
        v_ppfrag_HeHe*=             H_fraction;


       if(debug==1)
       {
        cout<<"process="<<process;
        if(process==1) cout<<" CR protons ";
        if(process==2) cout<<" CR He      ";
	cout<<" E_proton_total_rebin_GeV=" <<E_proton_total_rebin_GeV
            <<       " E_gamma_rebin_GeV=" <<E_gamma_rebin_GeV
            <<" comparing Huang="<<v<<" with Kachelriess2="<< v_ppfrag
            <<" K/H="<<v_ppfrag/v;
	if(process==1) cout<<" v_ppfrag_pHe=" <<v_ppfrag_pHe;
	if(process==2) cout<<" v_ppfrag_HeHe="<<v_ppfrag_HeHe;
        cout<<endl;

        if(process==1)cout<<" CR protons ";
        if(process==2)cout<<" CR He      ";
	cout<<" E_projectile =" <<E_proton_total_rebin_GeV
            <<" E_gamma=" <<E_gamma_rebin_GeV
            <<" comparing Huang="<<v<<" with Kachelriess2="<< v_ppfrag
            <<" K/H="<<v_ppfrag/v
            <<endl;
       }// if debug==1


        v= v_ppfrag;
	

     }// model_name=="Kachelriess2" || "Kamae_Kachelriess" || "Kamae"


       ////////////////////////////////////////////////////////////

       // Dermer, copy of Kachelriess2 with Dermer replacing Kamae

       if(model_name=="Dermer_Kachelriess" ||  model_name=="Dermer")
       {
	 // using combination of all  Dermer  with QCSJET
        double ep,es;
        int id=0;  // gammas
        int reac;  // reaction
        double ethr=20.;   // threshold in GeV used  for Dermer/QCSJET-II switch. Should be parameter. 
        if( model_name=="Dermer") ethr=1e6;// use Dermer only so set threshold arbitrarily high 
	

        ep = E_proton_total_rebin_GeV  ;
        es = E_gamma_rebin_GeV;

	double v_ppfrag,v_ppfrag_pHe, v_ppfrag_HeHe;;
        double A1,A2,A_factor;

        double  H_fraction=0.9;// H_fraction defined for Huang, in function emissivity, use the same here
        double He_fraction=1.0-H_fraction;
       
  


       if(process==1) // CR protons
       {

           
        
        if(ep<ethr) 
	{
	  // Dermer function uses kinetic energy in GeV
          // returns cross-section in mb/GeV times gamma-ray energy in GeV
	
         double Tp=ep-m_p; //       notation for kinetic energy of proton in GeV


  
        // use this notation for consistency with other  sections. 
         v_ppfrag=  EdsigDermerdE(es, Tp); //Esec, Tp in GeV, Units of mb

         if(debug==1)
         cout<<"Dermer routines: CR p  Tp="<<Tp<<" es="<<" v_ppfrag="<<v_ppfrag<<endl;


	 A1=1.;
         A2=4.;
	 A_factor=pow(A1*A2, 0.8); // Norbury, J. W., & Townsend, L. W. 2007, Nucl. Instrum. Meth. B 254, 187 
         v_ppfrag_pHe = v_ppfrag * A_factor * He_fraction;

	}

        if (ep>=ethr)
	{
         reac=1; // pp QCSJET-II, >10 GeV only

         v_ppfrag   = ppfrag_spec_int(ep,es,id,reac);

         reac=2; // p-He QCSJET-II >10 GeV only

         v_ppfrag_pHe   = ppfrag_spec_int(ep,es,id,reac);
         v_ppfrag_pHe  *= He_fraction;


        }

          v_ppfrag += v_ppfrag_pHe;
	
       }// process==1

  
  

	if(process==2) // CR He
      {


        if(ep<ethr)
	{
         double Tp=ep-m_p; //       notation for kinetic energy of proton in GeV

 
   

        // use this notation for consistency with rest of Kachelriess section
        v_ppfrag=  EdsigDermerdE(es, Tp); //Esec, Tp in GeV, Units of mb

	 if(debug==1)
         cout<<"Dermer routines: CR He: Tp="<<Tp<<" es="<<es<<" v_ppfrag="<<v_ppfrag<<endl;

         v_ppfrag_HeHe= v_ppfrag;

	 A1=4.;
         A2=1.;
	 A_factor=pow(A1*A2, 0.8); // Norbury, J. W., & Townsend, L. W. 2007, Nucl. Instrum. Meth. B 254, 187 
         v_ppfrag  *= A_factor; // approximate correction for He projectile

	 // approximation for HeHe
	 A1=4.;
         A2=4.;
	 A_factor=pow(A1*A2, 0.8);  
         v_ppfrag_HeHe *= A_factor *He_fraction;

         v_ppfrag += v_ppfrag_HeHe;
	}

        if(ep>=ethr)
	{
         reac=3; // He-p, >=10 GeV only
        
         v_ppfrag   = ppfrag_spec_int(ep,es,id,reac);

         v_ppfrag_HeHe=v_ppfrag;

	 // approximate correction for HeHe, relative to He-p
	 A1=4.;
         A2=4.;
	 A_factor=pow(A1*A2, 0.8)/pow(A1*1.0, 0.8);  
         v_ppfrag_HeHe *= A_factor * He_fraction;

         v_ppfrag += v_ppfrag_HeHe;

	}      

      } // process==2

     // returns dsigma/dlogE = E dsigma/dE, so divide by E. Also for other values for printed comparison

	v_ppfrag       /= es;
	v_ppfrag_pHe   /= es;
	v_ppfrag_HeHe  /= es;

        v_ppfrag     *=             H_fraction; // to correspond to Huang per atom definition assumed in function emissivity
        v_ppfrag_pHe *=             H_fraction;
        v_ppfrag_HeHe*=             H_fraction;


       if(debug==1)
       {	
        cout<<"process="<<process;
        if(process==1) cout<<" CR protons ";
        if(process==2) cout<<" CR He      ";
	cout<<" E_proton_total_rebin_GeV=" <<E_proton_total_rebin_GeV
            <<       " E_gamma_rebin_GeV=" <<E_gamma_rebin_GeV
            <<" comparing Huang="<<v<<" with Dermer_Kachelriess="<< v_ppfrag
            <<" D/H="<<v_ppfrag/v;
	if(process==1) cout<<" v_ppfrag_pHe=" <<v_ppfrag_pHe;
	if(process==2) cout<<" v_ppfrag_HeHe="<<v_ppfrag_HeHe;
        cout<<endl;

        if(process==1)cout<<" CR protons ";
        if(process==2)cout<<" CR He      ";
	cout<<" E_projectile =" <<E_proton_total_rebin_GeV
            <<" E_gamma=" <<E_gamma_rebin_GeV
            <<" comparing Huang="<<v<<" with Dermer_Kachelriess="<< v_ppfrag
            <<" D/H="<<v_ppfrag/v
            <<endl;
       }//if debug==1

        v= v_ppfrag;
	

     }// model_name=="Dermer_Kachelriess" || "Dermer"




      ////////////////////////////////////////////////////////////

       // Kachelriess and Ostapchenko ppfrag v02 2013

       if(model_name=="Kachelriess2013" )
       {
	 // using interpolation QGSJET/Kamae
        double ep,es;
        int id=0;  // gammas
        int reac;  // reaction
        int kamae_QGS = 1; // Kamae-QGSJET interpolation
	

        ep = E_proton_total_rebin_GeV  ;
        es = E_gamma_rebin_GeV;

	double v_ppfrag,v_ppfrag_pHe, v_ppfrag_HeHe;;
        double A1,A2,A_factor;

        double  H_fraction=0.9;// H_fraction defined for Huang, in function emissivity, use the same here
        double He_fraction=1.0-H_fraction;
       
     


       if(process==1) // CR protons
       {

           

	{
         reac=1; // pp

         v_ppfrag   = ppfrag_v02(ep,es,id,reac,kamae_QGS);

         reac=2; // p-He 

         v_ppfrag_pHe   = ppfrag_spec_int(ep,es,id,reac);
         v_ppfrag_pHe  *= He_fraction;


        }

          v_ppfrag += v_ppfrag_pHe;
	
       }// process==1

  
  

	if(process==2) // CR He
      {



	//       if(ep>=ethr)
	{
         reac=3; // He-p
        
         v_ppfrag   = ppfrag_v02(ep,es,id,reac,kamae_QGS);

         reac=4; // He-He 

         v_ppfrag_HeHe   =  ppfrag_v02(ep,es,id,reac,kamae_QGS);
         v_ppfrag_HeHe  *= He_fraction;

         v_ppfrag += v_ppfrag_HeHe;

	}      

      } // process==2

     // returns dsigma/dlogE = E dsigma/dE, so divide by E. Also for other values for printed comparison

	v_ppfrag       /= es;
	v_ppfrag_pHe   /= es;
	v_ppfrag_HeHe  /= es;

        v_ppfrag     *=             H_fraction; // to correspond to Huang per atom definition assumed in function emissivity
        v_ppfrag_pHe *=             H_fraction;
        v_ppfrag_HeHe*=             H_fraction;


       if(debug==1)
       {	
        cout<<"process="<<process;
        if(process==1) cout<<" CR protons ";
        if(process==2) cout<<" CR He      ";
	cout<<" E_proton_total_rebin_GeV=" <<E_proton_total_rebin_GeV
            <<       " E_gamma_rebin_GeV=" <<E_gamma_rebin_GeV
            <<" v_ppfrag="<< v_ppfrag;
            
	if(process==1) cout<<" v_ppfrag_pHe=" <<v_ppfrag_pHe;
	if(process==2) cout<<" v_ppfrag_HeHe="<<v_ppfrag_HeHe;
        cout<<endl;

 
       }//if debug==1

        v= v_ppfrag;
	

     }// model_name=="Kachelriess2013




       //////////////////////////////////////////////////////////////////////
    

       // Dermer combined with Kachelriess and Ostapchenko ppfrag v02 2013

       if(model_name=="Dermer_Kachelriess2013" || model_name=="Dermer_Kachelriess2013_heavy" )
       {
	 // using interpolation QGSJET/Kamae
        double ep,es;
        int id=0;  // gammas
        int reac;  // reaction
        int kamae_QGS = 1; // Kamae-QGSJET interpolation
	

        ep = E_proton_total_rebin_GeV  ;
        es = E_gamma_rebin_GeV;

	double v_ppfrag,v_ppfrag_pHe, v_ppfrag_HeHe;;
        double A1,A2,A_factor;
	double heavies_factor;

        double  H_fraction=0.9;// H_fraction defined for Huang, in function emissivity, use the same here
        double He_fraction=1.0-H_fraction;
       
        double ethr=20.;   // threshold in GeV used  for Dermer/QCSJET-II switch. Should be parameter. 

        int test_vary_heavies=1; // to see effect at low energy. normally this parameter should be 0

                           double A_factor_power=0.8 ;// Norbury, J. W., & Townsend, L. W. 2007, Nucl. Instrum. Meth. B 254, 187
        if( test_vary_heavies==1) A_factor_power=0.93;// extreme case: use this for all energies

       if(process==1) // CR protons
       {

           
       if(ep<=ethr) 
	{
	  // Dermer function uses kinetic energy in GeV
          // returns cross-section in mb/GeV times gamma-ray energy in GeV
	
         double Tp=ep-m_p; //       notation for kinetic energy of proton in GeV


  
        // use this notation for consistency with other  sections. 
         v_ppfrag=  EdsigDermerdE(es, Tp); //Esec, Tp in GeV, Units of mb

         if(debug==1)
         cout<<"Dermer routines: CR p  Tp="<<Tp<<" es="<<" v_ppfrag="<<v_ppfrag<<endl;


	 A1=1.;
         A2=4.;
	 A_factor=pow(A1*A2,  A_factor_power); // Norbury, J. W., & Townsend, L. W. 2007, Nucl. Instrum. Meth. B 254, 187 
         v_ppfrag_pHe = v_ppfrag * A_factor * He_fraction;


	   // HadronicHeaviesFactor.cc: power=0.8 returned factor (p+>He)*all  / (pp+pHe)  =1.06558
                  heavies_factor=1.06558;

		  if( test_vary_heavies==1) heavies_factor=1.09561;


	}//      ep<=ethr 

      if(ep>ethr) 
      {
         reac=1; // pp

         v_ppfrag   = ppfrag_v02(ep,es,id,reac,kamae_QGS);

         reac=2; // p-He 

         v_ppfrag_pHe   = ppfrag_spec_int(ep,es,id,reac);
         v_ppfrag_pHe  *= He_fraction;

	   // HadronicHeaviesFactor.cc: power=0.93 returned factor (p+>He)*all  / (pp+pHe)  =1.09561
                  heavies_factor==1.09561;


       } //     ep>ethr

          v_ppfrag += v_ppfrag_pHe;
	
	  if( model_name=="Dermer_Kachelriess2013_heavy")
	  {

	    //cout<<"process="<<process<<": using  heavies_factor="<<heavies_factor<<endl;
	   v_ppfrag*=heavies_factor;
	  }

       }// process==1

  
  

	if(process==2) // CR He
      {

      if(ep<=ethr)
	{
         double Tp=ep-m_p; //       notation for kinetic energy of proton in GeV

 
   

        // use this notation for consistency with rest of Kachelriess section
        v_ppfrag=  EdsigDermerdE(es, Tp); //Esec, Tp in GeV, Units of mb

	 if(debug==1)
         cout<<"Dermer routines: CR He: Tp="<<Tp<<" es="<<es<<" v_ppfrag="<<v_ppfrag<<endl;

         v_ppfrag_HeHe= v_ppfrag;

	 A1=4.;
         A2=1.;
	 A_factor=pow(A1*A2,  A_factor_power); // Norbury, J. W., & Townsend, L. W. 2007, Nucl. Instrum. Meth. B 254, 187 
         v_ppfrag  *= A_factor; // approximate correction for He projectile

	 // approximation for HeHe
	 A1=4.;
         A2=4.;
	 A_factor=pow(A1*A2,  A_factor_power);  
         v_ppfrag_HeHe *= A_factor *He_fraction;

         v_ppfrag += v_ppfrag_HeHe;

	   // power=0.8  HadronicHeaviesFactor.cc:returned factor     He *all  / (Hep+HeHe)=1.00736
                  heavies_factor=1.00736;


		  if( test_vary_heavies==1) heavies_factor=1.01026;


	}//ep<=ethr


	if(ep>=ethr)
	{
         reac=3; // He-p
        
         v_ppfrag   = ppfrag_v02(ep,es,id,reac,kamae_QGS);

         reac=4; // He-He 

         v_ppfrag_HeHe   =  ppfrag_v02(ep,es,id,reac,kamae_QGS);
         v_ppfrag_HeHe  *= He_fraction;

         v_ppfrag += v_ppfrag_HeHe;

	   // HadronicHeaviesFactor.cc: power=0.93 returned factor     He *all  / (Hep+HeHe)=1.01026
                  heavies_factor=1.01026;

	}//  (ep>ethr    

	  if( model_name=="Dermer_Kachelriess2013_heavy")
	  {

	    //cout<<"process="<<process<<": using  heavies_factor="<<heavies_factor<<endl;
	   v_ppfrag*=heavies_factor;
	  }


      } // process==2

     // returns dsigma/dlogE = E dsigma/dE, so divide by E. Also for other values for printed comparison

	v_ppfrag       /= es;
	v_ppfrag_pHe   /= es;
	v_ppfrag_HeHe  /= es;

        v_ppfrag     *=             H_fraction; // to correspond to Huang per atom definition assumed in function emissivity
        v_ppfrag_pHe *=             H_fraction;
        v_ppfrag_HeHe*=             H_fraction;


       if(debug==1)
       {	
        cout<<"process="<<process;
        if(process==1) cout<<" CR protons ";
        if(process==2) cout<<" CR He      ";
	cout<<" E_proton_total_rebin_GeV=" <<E_proton_total_rebin_GeV
            <<       " E_gamma_rebin_GeV=" <<E_gamma_rebin_GeV
            <<" v_ppfrag="<< v_ppfrag;
            
	if(process==1) cout<<" v_ppfrag_pHe=" <<v_ppfrag_pHe;
	if(process==2) cout<<" v_ppfrag_HeHe="<<v_ppfrag_HeHe;
        cout<<endl;

 
       }//if debug==1

        v= v_ppfrag;
	

     }// model_name=="Dermer_Kachelriess2013




       //////////////////////////////////////////////////////////////////////








       //////////////////////////////////////////////////////////////////////
       
  if(model_name=="FLUKA")
  {

     double E_proton_rebin_GeV=E_proton_rebin[iE_proton_rebin];
     if(energy_units=="MeV")  E_proton_rebin_GeV*= 1.e-3;         // MeV -> GeV

     iE_proton=(log(E_proton_rebin_GeV) - log(E_proton_FLUKA[0]))/log(E_proton_FLUKA[1]/ E_proton_FLUKA[0] ) +.0001; //  FLUKA uses kinetic energy and GeV
     iE_gamma =(log(E_gamma_rebin_GeV ) - log(E_gamma_FLUKA [0]))/log(E_gamma_FLUKA [1]/ E_gamma_FLUKA [0] ) +.0001; //  FLUKA uses                    GeV

     if(debug==1)
      {
      cout<<"FLUKA rebinning: ";
      cout<<"E_proton_rebin="<<E_proton_rebin[ iE_proton_rebin]<<" "<<energy_units<<" iE_proton="
		      << iE_proton<<" E_proton="<<E_proton_FLUKA[iE_proton]<<" - "<<E_proton_FLUKA[iE_proton+1] ;
     
      cout<<" E_gamma_rebin="<<E_gamma_rebin[ iE_gamma_rebin]<<" "<<energy_units<<" iE_gamma="
		      << iE_gamma<<" E_gamma="<<E_gamma_FLUKA[iE_gamma]<<" - "<<E_gamma_FLUKA[iE_gamma+1] << " GeV";
      }


     if(iE_proton>=0 && iE_proton<nE_proton_FLUKA-1)
     if(iE_gamma >=0 && iE_gamma <nE_gamma_FLUKA -1)
     {
	double v1= gamma_p_FLUKA_matrix[iE_proton  ][iE_gamma];
	double v2= gamma_p_FLUKA_matrix[iE_proton+1][iE_gamma];

	double v3= gamma_p_FLUKA_matrix[iE_proton  ][iE_gamma+1];
	double v4= gamma_p_FLUKA_matrix[iE_proton+1][iE_gamma+1];

 
	// first interpolate to E_proton_rebin for each E_gamma
        // differential cross-section (mb GeV-1)  

        double v12= v1+
                   (v2-v1)
                  *(E_proton_rebin_GeV                    - E_proton_FLUKA[ iE_proton])        
                  /(E_proton_FLUKA [iE_proton+1]          - E_proton_FLUKA[ iE_proton]);



        double v34= v3+
                   (v4-v3)
                  *(E_proton_rebin_GeV                    - E_proton_FLUKA[ iE_proton])         
                  /(E_proton_FLUKA[iE_proton+1]           - E_proton_FLUKA[ iE_proton]);

        
	if(debug==1) cout<<" v before reassignment="<<v;


	// then interpolate to E_gamma_rebin

               v=  v12+
	          (v34-v12)
                  *(E_gamma_rebin_GeV                   - E_gamma_FLUKA [ iE_gamma ])
                  /(E_gamma_FLUKA[iE_gamma +1]          - E_gamma_FLUKA [ iE_gamma ]);


	       if(debug==1) cout<<" v1="<<v1<<" v2="<<v2<<" v12="<<v12<<" v3="<<v3<<" v4="<<v4<<" v34="<<v34<<" v="<<v<<endl;

     } // if iE_proton iE_gamma


        double A1,A2,A_factor;

        double  H_fraction=0.9;// H_fraction defined for Huang, in function emissivity, use the same here
        double He_fraction=1.0-H_fraction;

	if(process==1)// CR protons + p + He
	{
	 A1=1.;
         A2=4.;
	 A_factor=pow(A1*A2, 0.8); // Norbury, J. W., & Townsend, L. W. 2007, Nucl. Instrum. Meth. B 254, 187 
         v += v*  A_factor * He_fraction;
	}


	if(process==2)// CR He + p
	{
	 A1=4.;
         A2=1.;
	 A_factor=pow(A1*A2, 0.8); // Norbury, J. W., & Townsend, L. W. 2007, Nucl. Instrum. Meth. B 254, 187 
         double vHep= v * A_factor;
         

	 // approximation for CR He + He
	 A1=4.;
         A2=4.;
	 A_factor=pow(A1*A2, 0.8);
	

         double vHeHe=v * A_factor*He_fraction;
         v =  vHep +  vHeHe;
      
	}

     }// model_name=="FLUKA"
       
       //////////////////////////////////////////////////////////////////////////////////////////////////////////////



        if(debug==1) cout<<"v before gamma energy conversion ="<<v<<endl;

        if (energy_units=="MeV") v*=1.e-3;// gammas GeV-1 -> MeV-1 

	//     if(v>0.) cout<<"ProductionMatrics::rebin: model_name="<<model_name<<" v after gamma energy conversion ="<<v<<endl;
        if(debug==1) cout<<"ProductionMatrics::rebin: model_name="<<model_name<<" v after gamma energy conversion ="<<v<<endl;

        matrix_rebin[iE_proton_rebin][iE_gamma_rebin] = v;





	//    }//if  iE_gamma moved

     }// iE_gamma_rebin

     //    }//if  iE_proton moved

   } // if total energy>mp

  }// iE_proton_rebin




    if(debug==2)
    {
     cout<<" rebin: debug of               matrix_rebin="<<endl;

     for(int iE_gamma_rebin =0;iE_gamma_rebin <nE_gamma_rebin ;iE_gamma_rebin ++)
     {
      cout<<"iE_gamma_rebin="<<iE_gamma_rebin<<" E_gamma_rebin="<<E_gamma_rebin[iE_gamma_rebin]<<endl;

      for(int iE_proton_rebin=0;iE_proton_rebin<nE_proton_rebin;iE_proton_rebin++)
         cout<<                matrix_rebin[iE_proton_rebin][iE_gamma_rebin]<<" ";
    
      cout<<endl;
     }
    } 

  if(process==1)gamma_decay_p_matrix_rebin = matrix_rebin;// NB works right for multiple calls only if clear() called
  if(process==2)gamma_decay_he_matrix_rebin= matrix_rebin;// NB works right for multiple calls only if clear() called

  }// process


 if(debug==2)
    {
  cout<<" rebin: debug of gamma_decay_p_matrix_rebin="<<endl;

  for(int iE_gamma_rebin =0;iE_gamma_rebin <nE_gamma_rebin ;iE_gamma_rebin ++)
  {
    cout<<"iE_gamma_rebin="<<iE_gamma_rebin<<" E_gamma_rebin="<<E_gamma_rebin[iE_gamma_rebin]<<endl;

  for(int iE_proton_rebin=0;iE_proton_rebin<nE_proton_rebin;iE_proton_rebin++)
         cout<<  gamma_decay_p_matrix_rebin[iE_proton_rebin][iE_gamma_rebin]<<" ";
    
  cout<<endl;
  }

  cout<<" rebin: debug of gamma_decay_he_matrix_rebin="<<endl;

  for(int iE_gamma_rebin =0;iE_gamma_rebin <nE_gamma_rebin ;iE_gamma_rebin ++)
  {
    cout<<"iE_gamma_rebin="<<iE_gamma_rebin<<" E_gamma_rebin="<<E_gamma_rebin[iE_gamma_rebin]<<endl;

  for(int iE_proton_rebin=0;iE_proton_rebin<nE_proton_rebin;iE_proton_rebin++)
         cout<<  gamma_decay_he_matrix_rebin[iE_proton_rebin][iE_gamma_rebin]<<" ";
    
  cout<<endl;
  }
    }


  matrices_rebinned=1;

 cout<<"<<ProductionMatrices::rebin"<<endl;

  return 0;
};

///////////////////////////////////////////////////////////////////

  int ProductionMatrices::emissivity(valarray<double> proton_spectrum_rebin, valarray<double> helium_spectrum_rebin )
{
  if(debug==1) cout<<"ProductionMatrices::emissivity    "<<endl;

  // proton and helium spectra are assumed to be fluxes cm-2 sr-1 s-1 GeV-1
  // (so velocity term is already included)
  // giving emissivity sr-1 s-1 GeV-1 H atom-1

  if(matrices_rebinned!= 1) {cout<<" ProductionMatrices:: emissivity: matrices have not been not rebinned!"<<endl; exit(1);}

  gamma_decay_p_emissivity_rebin .resize(nE_gamma_rebin);//AWS20120328
  gamma_decay_he_emissivity_rebin.resize(nE_gamma_rebin);//AWS20120328

  gamma_decay_p_emissivity_rebin =0.;
  gamma_decay_he_emissivity_rebin=0.;

  int debug=0;

  double m_p   = 0.938 ;      // proton mass in GeV  (deprecated hard coding of physics constant !)
  double c     = 3.0e10;      // velocity of light   (deprecated hard coding  of physics constant!)
  
  // it is assumed the spectrum is on a log energy scale, hence multiply by E_proton and apply factor 

for (int iE_gamma_rebin =0; iE_gamma_rebin <nE_gamma_rebin ; iE_gamma_rebin++ )
  {
   for(int iE_proton_rebin=0; iE_proton_rebin<nE_proton_rebin; iE_proton_rebin++)
   {
     double g    = E_proton_rebin[iE_proton_rebin]/m_p;// Lorentz factor of proton
     double beta = sqrt(1.0-1.0/g/g);                  // velocity/c     of proton

     // energy is total energy per nucleon in GeV
     // since energy is per nucleon, protons and helium have same velocity
     // but velocity not needed if spectra are already fluxes. Keep for case of CR densities e.g. for GALPROP

    gamma_decay_p_emissivity_rebin[ iE_gamma_rebin]+= 
                                                        proton_spectrum_rebin     [iE_proton_rebin]
                                                      * gamma_decay_p_matrix_rebin[iE_proton_rebin][iE_gamma_rebin]
                                                      * E_proton_rebin[iE_proton_rebin];

    gamma_decay_he_emissivity_rebin[ iE_gamma_rebin]+=  
                                                        helium_spectrum_rebin      [iE_proton_rebin]
                                                      * gamma_decay_he_matrix_rebin[iE_proton_rebin][iE_gamma_rebin]
                                                      * E_proton_rebin[iE_proton_rebin];

 
    if(debug==1)
      cout<<"emissivity::test iE_gamma_rebin="<<iE_gamma_rebin  <<" E_gamma_rebin= "<< E_gamma_rebin [iE_gamma_rebin]
        <<" E_proton_rebin= "   << E_proton_rebin[iE_proton_rebin]
        <<" proton spectrum="<<proton_spectrum_rebin     [iE_proton_rebin]
	<<" matrix="<<gamma_decay_p_matrix_rebin[iE_proton_rebin][iE_gamma_rebin]
        <<" gamma_decay_p_emissivity_rebin= "<<gamma_decay_p_emissivity_rebin[ iE_gamma_rebin]
	<<" energy_units="<<energy_units
        <<endl;

  

    }
    }

  gamma_decay_p_emissivity_rebin  *=  log(E_proton_rebin[1]/E_proton_rebin[0]); // since log scale integration
  gamma_decay_he_emissivity_rebin *=  log(E_proton_rebin[1]/E_proton_rebin[0]); // since log scale integration



  double integral_proton_flux=0.;
  double integral_helium_flux=0.;

  for(int iE_proton_rebin=0; iE_proton_rebin<nE_proton_rebin; iE_proton_rebin++)
  {
   integral_proton_flux                          += 
                                                        proton_spectrum_rebin    [iE_proton_rebin]
                                                      * E_proton_rebin[iE_proton_rebin];

    integral_helium_flux                          += 
                                                        helium_spectrum_rebin    [iE_proton_rebin]
                                                      * E_proton_rebin[iE_proton_rebin];
}


  integral_proton_flux  *=  log(E_proton_rebin[1]/E_proton_rebin[0]);
  integral_helium_flux  *=  log(E_proton_rebin[1]/E_proton_rebin[0]);

  if(debug==1) cout<<"factor for log integration="<< log(E_proton_rebin[1]/E_proton_rebin[0])<<endl;

  gamma_decay_p_emissivity_rebin  *= 1.0e-27; // cross-section is in millibarn GeV-1
  gamma_decay_he_emissivity_rebin *= 1.0e-27; // cross-section is in millibarn GeV-1


  double H_fraction=0.9;  // Huang etal is per atom of ISM with H fraction 0.9. Convert to per H atom standard
  gamma_decay_p_emissivity_rebin /= H_fraction;
  gamma_decay_he_emissivity_rebin/= H_fraction;

  if(debug==1) cout<<"integral proton flux="<<integral_proton_flux<<" cm-2 sr-1 s-1"<<endl;
  if(debug==1) cout<<"integral helium flux="<<integral_helium_flux<<" cm-2 sr-1 s-1"<<endl;

  return 0;
}

///////////////////////////////////////////////////////////////////

valarray<double> ProductionMatrices::emissivity_total
                  (valarray<double> proton_spectrum_rebin, valarray<double> helium_spectrum_rebin )
{
  if(debug==1) cout<<"ProductionMatrices::emissivity_total    "<<endl;
  
    
 emissivity( proton_spectrum_rebin, helium_spectrum_rebin);
 return  gamma_decay_p_emissivity_rebin +  gamma_decay_he_emissivity_rebin;

}
///////////////////////////////////////////////////////////////////

int ProductionMatrices::rebin_with_energy_dispersion(valarray<double> E_proton_rebin_, valarray<double> E_gamma_rebin_, valarray<double> E_gamma_true_)
{
  // computes the matrix accounting for Fermi-LAT energy dispersion
  // uses E_gamma_true to compute response to wide enough energy range, then energy-dispersed response on E_gamma_rebin

  cout<<">> ProductionMatrices:::rebin_with_energy_dispersion"<<endl;

  cout<<"  ProductionMatrices:::rebin_with_energy_dispersion E_gamma_rebin_="; for(int i=0;i<E_gamma_rebin_.size();i++) cout<<E_gamma_rebin_[i]<<" "; cout<<endl;
  cout<<"  ProductionMatrices:::rebin_with_energy_dispersion E_gamma_true_ ="; for(int i=0;i<E_gamma_true_ .size();i++) cout<<E_gamma_true_ [i]<<" "; cout<<endl;


  // extended range matrix
  rebin( E_proton_rebin_, E_gamma_true_);

  // cout<<endl<<"  ProductionMatrices:::rebin_with_energy_dispersion: matrix for E_gamma_true :"<<endl;  print();
 
  cout<<"  ProductionMatrices:::rebin_with_energy_dispersion: energy dispersion starting"<<endl;
 int debug=0;

  int nE_meas=E_gamma_rebin_.size();
  valarray<double> work;

   // check integrals
  valarray<double> integral_true;
  valarray<double> integral_meas;
  integral_true.resize(nE_proton_rebin);
  integral_meas.resize(nE_proton_rebin);

  work.resize(E_gamma_true_.size());
  for (int iE_proton_rebin=0; iE_proton_rebin<nE_proton_rebin; iE_proton_rebin++)
  {
   work =  gamma_decay_p_matrix_rebin [iE_proton_rebin] * E_gamma_true_;
   integral_true[iE_proton_rebin]=work.sum()*log(E_gamma_true_[1]/E_gamma_true_[0]);

   if(debug==1)
   {
   cout<<endl<<"  ProductionMatrices:::rebin_with_energy_dispersion: matrix for E_gamma_true :";
   for(int i=0;i<E_gamma_true_.size();i++)   cout<<gamma_decay_p_matrix_rebin [iE_proton_rebin][i]<<" ";   cout<<endl;
   cout<<" ProductionMatrices:::rebin_with_energy_dispersion: iE_proton_rebin="<<iE_proton_rebin<<" integral true= "<<integral_true[iE_proton_rebin]  <<endl;
   }


  }



  work.resize(nE_meas);
 
  

  for (int iE_proton_rebin=0; iE_proton_rebin<nE_proton_rebin; iE_proton_rebin++)
  {
   energyDispersion.ApplyEnergyDispersion( E_gamma_true_, gamma_decay_p_matrix_rebin [iE_proton_rebin] ,  E_gamma_rebin_,  work,  debug);

   gamma_decay_p_matrix_rebin [iE_proton_rebin].resize(nE_meas);
   gamma_decay_p_matrix_rebin [iE_proton_rebin] = work;

   energyDispersion.ApplyEnergyDispersion( E_gamma_true_, gamma_decay_he_matrix_rebin[iE_proton_rebin] ,  E_gamma_rebin_,  work,  debug);

   gamma_decay_he_matrix_rebin[iE_proton_rebin].resize(nE_meas);
   gamma_decay_he_matrix_rebin[iE_proton_rebin] = work;

  }

  cout<<"  ProductionMatrices:::rebin_with_energy_dispersion: energy dispersion completed"<<endl;


  // reset the gamma grid  E_gamma_rebin in the class object to correspond to the modified grid
  nE_gamma_rebin =  nE_meas;
   E_gamma_rebin.resize(nE_gamma_rebin);
   E_gamma_rebin = E_gamma_rebin_;
   cout<<" ProductionMatrices::rebin_with_energy_dispersion: nE_meas="<<nE_meas<<" nE_gamma_rebin="<<nE_gamma_rebin<<endl;


   // check integrals
   work.resize(E_gamma_rebin_.size());
  for (int iE_proton_rebin=0; iE_proton_rebin<nE_proton_rebin; iE_proton_rebin++)
  {
   work =  gamma_decay_p_matrix_rebin [iE_proton_rebin] * E_gamma_rebin_;
   integral_meas[iE_proton_rebin]=work.sum()*log(E_gamma_rebin_[1]/E_gamma_rebin_[0]);

   if(debug==1)
   {
   cout<<" ProductionMatrices::rebin_with_energy_dispersion: matrix meas=";
   for(int i=0;i<nE_gamma_rebin;i++)
   cout<<gamma_decay_p_matrix_rebin [iE_proton_rebin][i]<<" ";
   cout<<endl;
   }

   cout<<" ProductionMatrices:::rebin_with_energy_dispersion: iE_proton_rebin="<<iE_proton_rebin
       <<" integral true= "<<integral_true[iE_proton_rebin]
       <<" meas="          <<integral_meas[iE_proton_rebin]
       <<" meas/true="     <<integral_meas[iE_proton_rebin]/integral_true[iE_proton_rebin]
       <<endl;
  }


  cout<<endl<<"  ProductionMatrices:::rebin_with_energy_dispersion: matrix for E_gamma_rebin_ :"<<endl;
  //  print();


  cout<<"<< ProductionMatrices:::rebin_with_energy_dispersion"<<endl;
  return 0;
}

///////////////////////////////////////////////////////////////////

double  ProductionMatrices::emissivity_band
          (valarray<double> proton_spectrum_rebin, valarray<double> helium_spectrum_rebin,
           valarray<double> E_proton_rebin_,
           double E_gamma_band_min, double E_gamma_band_max, double E_gamma_factor)
{

  //  int debug=1; now in class set by init

  if(debug==1)cout<<"ProductionMatrices::emissivity_band"<<endl;

  valarray<double> E_gamma_rebin_work;
  

  int nE_gamma_rebin_=log(E_gamma_band_max/E_gamma_band_min)/log(E_gamma_factor)+1.001; // redone in rebin ? check

  E_gamma_rebin_work.resize(nE_gamma_rebin_);

  for(int iE_gamma_rebin=0;iE_gamma_rebin<nE_gamma_rebin_;iE_gamma_rebin++)
        E_gamma_rebin_work[iE_gamma_rebin] = E_gamma_band_min * pow(E_gamma_factor,iE_gamma_rebin);

  

  rebin(E_proton_rebin_,E_gamma_rebin_work);

  emissivity(proton_spectrum_rebin, helium_spectrum_rebin);

  double emiss_band=0;//zero'd AWS20120327

  for(int iE_gamma_rebin=0;iE_gamma_rebin<nE_gamma_rebin;iE_gamma_rebin++)
    {

      emiss_band+= gamma_decay_p_emissivity_rebin  [iE_gamma_rebin]*E_gamma_rebin[iE_gamma_rebin];
      emiss_band+= gamma_decay_he_emissivity_rebin [iE_gamma_rebin]*E_gamma_rebin[iE_gamma_rebin];

      if(debug==1)cout<<"emissivity_band: E_gamma= "<<E_gamma_rebin[iE_gamma_rebin]<<" "<<energy_units<<" "
                                                   << gamma_decay_p_emissivity_rebin [iE_gamma_rebin]<<" "
                                                   << gamma_decay_he_emissivity_rebin[iE_gamma_rebin]<<endl;
    }

  emiss_band*=log(E_gamma_rebin[1]/E_gamma_rebin[0]);

  if(debug==1)cout<<"emissivity band: emiss_band="<<emiss_band<<endl;

  return emiss_band;
};

///////////////////////////////////////////////////////////////////

valarray<double>  ProductionMatrices::emissivity_band
          (valarray<double> proton_spectrum_rebin, valarray<double> helium_spectrum_rebin,
           valarray<double> E_proton_rebin_,
           valarray<double>E_gamma_band_min, valarray<double> E_gamma_band_max, double E_gamma_factor)
{

  cout<<"ProductionMatrices::valarray<double> emissivity_band "<<endl;

  int nE_gamma_band=E_gamma_band_min.size();
  valarray<double> emiss_band;
  emiss_band.resize(nE_gamma_band);
  emiss_band=0;//zero'd AWS20120327
  
  for (int iE_gamma_band=0; iE_gamma_band< nE_gamma_band; iE_gamma_band++)
  {
  
   emiss_band[iE_gamma_band] = emissivity_band
                              ( proton_spectrum_rebin,helium_spectrum_rebin,
                                E_proton_rebin_,
                                E_gamma_band_min[iE_gamma_band], E_gamma_band_max[iE_gamma_band],E_gamma_factor);

  if(debug==1)
  cout<<"emissivity band valarray<double>:"
             <<" E_gamma_band_min="<< E_gamma_band_min[iE_gamma_band]
             <<" E_gamma_band_max="<< E_gamma_band_max[iE_gamma_band]
             <<" emiss_band="<<emiss_band[iE_gamma_band]<<endl;

  }//  iE_gamma_band

  return emiss_band;
};
//////////////////////////////////////////////////////////////////

 valarray<double> ProductionMatrices::  emissivity_band
                   (valarray<double> proton_spectrum_rebin, valarray<double> helium_spectrum_rebin)
 {  
 // version using gamma_band matrices precomputed with gen_gamma_band_matrices
 
  int debug=0; //AWS20120621

  if(debug==1) cout<<"ProductionMatrices::valarray<double> emissivity_band using precomputed matrices for bands "<<endl;//AWS20120621



  // matrix precomputed
  int nE_gamma_band= gamma_band_p_matrix [0].size();

  valarray<double> emiss_band;
  emiss_band.resize(nE_gamma_band);
  emiss_band=0;


 for (int iE_gamma_band=0; iE_gamma_band< nE_gamma_band; iE_gamma_band++)
  {

  for(int iE_proton_rebin=0; iE_proton_rebin<nE_proton_rebin; iE_proton_rebin++)
   {
     //double g    = E_proton_rebin[iE_proton_rebin]/m_p;// Lorentz factor of proton
     //double beta = sqrt(1.0-1.0/g/g);                  // velocity/c     of proton

     // energy is total energy per nucleon in GeV
     // since energy is per nucleon, protons and helium have same velocity
     // but velocity not needed if spectra are already fluxes. Keep for case of CR densities e.g. for GALPROP

    emiss_band[ iE_gamma_band ]+= 
                                                        proton_spectrum_rebin     [iE_proton_rebin]
                                                      * gamma_band_p_matrix  [iE_proton_rebin][iE_gamma_band]
                                                      * E_proton_rebin[iE_proton_rebin];

    emiss_band[ iE_gamma_band ]+=  
                                                        helium_spectrum_rebin      [iE_proton_rebin]
                                                      * gamma_band_he_matrix [iE_proton_rebin][iE_gamma_band]
                                                      * E_proton_rebin[iE_proton_rebin];

 
    if(debug>=2)
      cout<<"emissivity_band ::using precomputed matrices iE_gamma_band ="<<iE_gamma_band  
        <<" E_proton_rebin= "   << E_proton_rebin[iE_proton_rebin]
        <<" proton spectrum="<<proton_spectrum_rebin[iE_proton_rebin]
        <<" Helium spectrum="<<helium_spectrum_rebin[iE_proton_rebin]
	<<" matrix="<<gamma_band_p_matrix  [iE_proton_rebin][iE_gamma_band]
	<<" emiss_band= "<<emiss_band[ iE_gamma_band ]
	<<" energy_units="<<energy_units
        <<endl;

  

   }//iE_proton_rebin
  }// iE_gamma_band

  emiss_band *=  log(E_proton_rebin[1]/E_proton_rebin[0]); // since log scale integration
  
  emiss_band *= 1.0e-27; // cross-section is in millibarn GeV-1
  
  double H_fraction=0.9;  // Huang etal is per atom of ISM with H fraction 0.9. Convert to per H atom standard
  emiss_band /= H_fraction;
  
    if(debug>=1)
    for (int iE_gamma_band=0; iE_gamma_band< nE_gamma_band; iE_gamma_band++)
     cout<<"emissivity_band ::result using precomputed matrices iE_gamma_band ="<<iE_gamma_band  
	<<" emiss_band= "<<emiss_band[ iE_gamma_band ]
	<<" energy_units="<<energy_units
        <<endl;

  return emiss_band;
}
///////////////////////////////////////////////////////////////////

int ProductionMatrices::gen_gamma_band_matrices
               ( valarray<double> E_proton_rebin_,
                 valarray<double>E_gamma_band_min, valarray<double> E_gamma_band_max, double E_gamma_factor)
{

  cout<<"ProductionMatrices::gen_gamma_band_matrices "<<endl;



  //  int debug=1; now in class set by init

 int nE_gamma_band=E_gamma_band_min.size();

  gamma_band_p_matrix .resize(E_proton_rebin_.size());
  gamma_band_he_matrix.resize(E_proton_rebin_.size());

 //for(int iE_proton_rebin=0;iE_proton_rebin<nE_proton_rebin;iE_proton_rebin++) AWS20120405 nE_proton_rebin not yet defined

 for(int iE_proton_rebin=0;iE_proton_rebin< E_proton_rebin_.size() ;iE_proton_rebin++)
  {
   gamma_band_p_matrix [iE_proton_rebin].resize(nE_gamma_band);
   gamma_band_he_matrix[iE_proton_rebin].resize(nE_gamma_band);
  }

  
  if(debug>=1)
  for(int iE_gamma_band=0; iE_gamma_band< nE_gamma_band; iE_gamma_band++)
         cout<<"gen_gamma_band_matrices: input"
             <<" E_gamma_band_min="<< E_gamma_band_min[iE_gamma_band]
             <<" E_gamma_band_max="<< E_gamma_band_max[iE_gamma_band]
             <<endl;

  valarray<double> E_gamma_rebin_work;


 for(int iE_gamma_band=0; iE_gamma_band< nE_gamma_band; iE_gamma_band++)
 {

  int nE_gamma_rebin_=
   log(E_gamma_band_max[iE_gamma_band]/E_gamma_band_min[iE_gamma_band])/log(E_gamma_factor)+1.001;//redone in rebin ? check

  if(debug>=1) cout<<"ProductionMatrices::gen_gamma_band_matrices: nE_gamma_rebin_="<<nE_gamma_rebin_<<endl;

  E_gamma_rebin_work.resize(nE_gamma_rebin_);

  if(debug>=1) cout<<"ProductionMatrices::gen_gamma_band_matrices: after E_gamma_rebin_work.resize"<<endl;

  for(int iE_gamma_rebin=0;iE_gamma_rebin<nE_gamma_rebin_;iE_gamma_rebin++)
        E_gamma_rebin_work[iE_gamma_rebin] = E_gamma_band_min[iE_gamma_band] * pow(E_gamma_factor,iE_gamma_rebin);

  
  

  if (energy_dispersion==0)
  rebin(E_proton_rebin_,E_gamma_rebin_work); // NB resets E_proton_rebin, E_gamma_rebin in class.

  
  // test energy dispersion

  if (energy_dispersion==1)
  {
   double fac_min=0.5;
   double fac_max=2;
   double E_gamma_factor_true = 1+(E_gamma_factor-1)*0.5;
   //       fac_min=1.0;fac_max=1.0;// to test
   //       fac_min=1.2;fac_max=0.8;// to test
   double E_gamma_true_min= E_gamma_band_min [iE_gamma_band] * fac_min;
   double E_gamma_true_max= E_gamma_band_max [iE_gamma_band] * fac_max;
 
 
   int nE_gamma_true=  log(E_gamma_true_max/E_gamma_true_min)/log(E_gamma_factor_true)+1.001;
   valarray<double> E_gamma_true;
   E_gamma_true.resize( nE_gamma_true);
                
                  
   for(int i=0;i<nE_gamma_true;i++) E_gamma_true[i]=   E_gamma_true_min   *  pow(E_gamma_factor_true ,i);
   rebin_with_energy_dispersion(E_proton_rebin_, E_gamma_rebin_work, E_gamma_true  );
  }


  // from here unchanged

  for(int iE_proton_rebin=0;iE_proton_rebin<nE_proton_rebin;iE_proton_rebin++)
  {
   double emiss_band_p=0;
   double emiss_band_he=0;
  
   for(int iE_gamma_rebin=0;iE_gamma_rebin<nE_gamma_rebin;iE_gamma_rebin++)
   {

      emiss_band_p += gamma_decay_p_matrix_rebin [iE_proton_rebin] [iE_gamma_rebin]*E_gamma_rebin[iE_gamma_rebin];
      emiss_band_he+= gamma_decay_he_matrix_rebin[iE_proton_rebin] [iE_gamma_rebin]*E_gamma_rebin[iE_gamma_rebin];
      

      if(debug>=3)cout<<"gen_gamma_band_matrices: "<<" E_proton_rebin="<<E_proton_rebin [iE_proton_rebin]
                                                   <<" E_gamma= "      <<E_gamma_rebin[iE_gamma_rebin]
                                                   <<" "<<energy_units
	                      	                   <<" emiss_band_p=" << emiss_band_p
		                                   <<" emiss_band_he="<< emiss_band_he
                                                   <<endl;
    }

  emiss_band_p *=log(E_gamma_rebin[1]/E_gamma_rebin[0]);
  emiss_band_he*=log(E_gamma_rebin[1]/E_gamma_rebin[0]);

 


	 gamma_band_p_matrix  [iE_proton_rebin][iE_gamma_band] = emiss_band_p;
	 gamma_band_he_matrix [iE_proton_rebin][iE_gamma_band] = emiss_band_he;

       if(debug>=2)
         cout<<"gen_gamma_band_matrices: result"
             <<" E_proton_rebin="  << E_proton_rebin  [iE_proton_rebin]
             <<" E_gamma_band_min="<< E_gamma_band_min[iE_gamma_band]
             <<" E_gamma_band_max="<< E_gamma_band_max[iE_gamma_band]
             <<" "<<energy_units
             <<" emiss_band_p="<<emiss_band_p
             <<" emiss_band_he="<<emiss_band_he
             <<endl;


    }// iE_gamma_band

  }// iE_proton_rebin





 if(debug>=1)
 for(int iE_proton_rebin=0;iE_proton_rebin<nE_proton_rebin;iE_proton_rebin++)
 for(int iE_gamma_band=0; iE_gamma_band< nE_gamma_band; iE_gamma_band++)

         if(debug==1)
         cout<<"gen_gamma_band_matrices: result matrices"
             <<" E_proton_rebin="  << E_proton_rebin  [iE_proton_rebin]
             <<" E_gamma_band_min="<< E_gamma_band_min[iE_gamma_band]
             <<" E_gamma_band_max="<< E_gamma_band_max[iE_gamma_band]
             <<" "<<energy_units
             <<" p: "<< gamma_band_p_matrix  [iE_proton_rebin][iE_gamma_band]
             <<" He:"<< gamma_band_he_matrix [iE_proton_rebin][iE_gamma_band]
             <<endl;



  return 0;
};


/////////////////////////////////////////////////////////////////
// Kachelriess & Ostapchenko 2012: ppfrag
////////////////////////////////////////////////////////////////


// g++ frag.f ppfrag.cc -lgfortran -lgfortranbegin

// for normal usage:
void   ppfrag_spec_ini();
double ppfrag_spec_int (double ,double ,int ,int );


///////////////////////////////////////////////////////////////
// interface to ppfrag fortran routines in frag.f
extern "C" void   spec_ini_();
extern "C" double spec_int_(double* ,double* ,int* ,int* );

// ppv02 AWS20130204
extern "C" void     init_();                                     // init.f95 : reads gam_pp gam_pHe gam_Hep gam_HeHe datasets
extern "C" void     spectrum_(double*, double*, int*, double*);  // in photons.f95:  subroutine spectrum(E_p,E_g,reac,spec)
extern "C" double   spec_qgs_(double*, double*, int*);           // in photons.f95:  double precision funcion  spec_QGS(E_p,E_g,k)   k=reac. NB lower case required !
///////////////////////////////////////////////////////////////

void   ProductionMatrices::ppfrag_spec_ini()
{

  // reads data from gamfrag.dat and apfrag.dat
  cout<<">> ProductionMatrices::ppfrag_spec_ini"<<endl;
  cout<<"   reading Kachelriess ppfrag matrices"<<endl;
  spec_ini_();

  cout<<"   ppv02 init_()"<<endl; // AWS20130204
  init_();

  cout<<"<< ProductionMatrices::ppfrag_spec_ini"<<endl;
  return;
}
///////////////////////////////////////////////////////////////
double ProductionMatrices::ppfrag_spec_int(double ep,double es,int id, int reac)
{
  /*
C++ interface to fortran routines described in
http://arxiv.org/abs/1206.4705 M. Kachelriess, S. Ostapchenko
http://sourceforge.net/projects/ppfrag/ Version 2012-06-21

Original header from frag.f:
!-----------------------------------------------------------------------------
! interpolation of lab. energy spectra of photons and antiprotons (+antineutrons)
! input:
!   ep  - incident energy per nucleon (GeV);
!   es  - secondary particle energy (GeV);
!   id  - secondary particle type (0 - photon, 1 - antiproton)
!   reac - production process (1 - p+p; 2 - p+He; 3 - He+p)
!   reac=0: photon spectra for p+p collisions based on the non-diffractive part
!       of the Kamae parametrization for ep<ethr and present tables for ep>ethr 
!-----------------------------------------------------------------------------


NB incident particle energy ep is TOTAL per nucleon (see frag.f)
returns differential cross-section in mb/GeV, multiplied by es in GeV (dsigma/dlogE)
From the cited paper, ethr=50 GeV. 

  */

  if(ppfrag_read!=1) { ppfrag_spec_ini(); ppfrag_read=1;}

  double result=0.;

  // respect limits and avoid stop in fortran routines

  // pp including Kamae
  if(reac==0 && es<ep && id==0)   result=spec_int_(&ep,&es,&id,&reac);



  // other cases including pbar, only >= 10 GeV
  if(reac >0 && es<ep && ep>=10.) result=spec_int_(&ep,&es,&id,&reac);


  int ppfragv02 =0; // used to compare versions. normally 0

  if( ppfragv02==1 && id==0)
  {
  //ppfrag v02
  // NB only defined for reac=1-4  1=p+p 2=p+He 3=He+p 4=He+He
  {
  if(id==0)cout<<" ProductionMatrices::ppfrag_spec_int: id="<<id<<" reac="<<reac<<" ep="<<ep<<" es="<<es<<" GeV, result from spec_int="<<result <<" mb/GeV*es; " ;

  double ep_eV=ep*1e9; // requires eV
  double es_eV=es*1e9; // requires eV

  if(reac==0) reac=1; // since this routine interpolates to Kamae, corresponding to original definition of reac=0
  spectrum_(&ep_eV,&es_eV,    &reac,&result); // result in mb/MeV (not mb/eV as written in photons.f95)
  result*=1.e3*es; // for compatibility with spec_int:   mb/GeV * es
  cout<<" result from ppfrag v02 spectrum times 1e3*es="<<result;

  result=spec_qgs_(&ep_eV,&es_eV,    &reac); // result in mb/GeV * es
  cout<<" result from ppfrag v02 spec_QGS ="<<result;
  cout<<endl;
  }
  }// ppfragv02==1

  return result;
}

///////////////////////////////////////////////////////////////////
double ProductionMatrices::ppfrag_v02(double ep,double es,int id, int reac, int kamae_QGS)
{
  /*
C++ interface to ppfrag v02 init.f95 photons.f95 routines from  M. Kachelriess, S. Ostapchenko
For earlier version see above and
http://arxiv.org/abs/1206.4705 M. Kachelriess, S. Ostapchenko
http://sourceforge.net/projects/ppfrag/ Version 2012-06-21

Interpolation to Kamae cross-sections at low energies is implemented here.

!-----------------------------------------------------------------------------
! 
! input:
!   ep  - incident energy per nucleon (GeV);
!   es  - secondary particle energy (GeV);
!  
!   reac - reaction:  1 = p+p; 2 = p+He; 3 = He+p; 4 = He+He 
!                                                                                                                                                               
!-----------------------------------------------------------------------------


NB incident particle energy ep is TOTAL per NUCLEON (for consistency with ppfrag_spec_int)
   as in frag.f which is different from photons.f95 (total per NUCLEUS) so have to convert for He
   returns differential cross-section in mb/GeV, multiplied by es in GeV (dsigma/dlogE)


  */

  int debug=0;

  if(ppfrag_read!=1) { ppfrag_spec_ini(); ppfrag_read=1;} // reads data for both versions

  double result=0.;
  
  //ppfrag v02
  // NB only defined for reac=1-4  1=p+p 2=p+He 3=He+p 4=He+He
  
  if(reac<1 || reac>4){cout<<" ProductionMatrices::ppfrag_v02 invalid reac ! :"<<reac<<endl;  return result;}

  //if(id==0)cout<<" ProductionMatrices::ppfrag_spec_int: id="<<id<<" reac="<<reac<<" ep="<<ep<<" es="<<es<<" GeV, result from spec_int="<<result <<" mb/GeV*es; " ;

  double ep_eV=ep*1e9; // requires eV
  double es_eV=es*1e9; // requires eV

  if(reac==3||reac==4) ep_eV*=4.0; // spectrum_ requires total energy per NUCLEUS

 if(kamae_QGS==1) // interpolated with Kamae
 {
  spectrum_(&ep_eV,&es_eV,    &reac,&result); // result in mb/MeV (not mb/eV as written in photons.f95)
  result*=1.e3*es; // for compatibility with spec_int 
 
  if(debug==1) cout<<" ProductionMatrices::ppfrag_v02: id="<<id<<" reac="<<reac<<" ep="<<ep<<" es="<<es<<" GeV, result="<<result <<" mb/GeV*es"<<endl; 
 }

 if(kamae_QGS==0) // QGS only
 {
  result=spec_qgs_(&ep_eV,&es_eV,    &reac); // result in mb/GeV * es
  cout<<" result from ppfrag_v02 spec_QGS ="<<result<<endl;
 }




  if(isnan(result) || isinf(result))
  {
   cout<<" ProductionMatrices::ppfrag_v02: id="<<id<<" reac="<<reac<<" ep="<<ep<<" es="<<es<<" GeV: result is nan or inf, setting to 0 !"<<endl;
   result=0.; // since returns nan for too large ep and  too small es  
  }

  if(isnan(result)||isinf(result)) exit(-1); // to make sure

  return result;
}
////////////////////////////////////////////////////

int   ProductionMatrices::test_ppfrag()
{
  cout<<"ppfrag.cc test program"<<endl;

  //  ppfrag_spec_ini(); done authomatically if required

  double ep,es;
  int id,reac;
  int kamae_QGS;// 0=QGSJET only 1=interpolation

  ep=1000.;// proton energy GeV
  es=   1.;// gamma  energy GeV

  double epfac =pow(10,0.1);// 10 per decade
  double esfac =pow(10,0.1);// 10 per decade

  id=0;  // 0=photons 1=antiprotons
  reac=0;// reaction

  double result0,result1,   result2,   result3;
  double         result1v02,result2v02,result3v02,result4v02; //AWS20130205

  cout<<"============ ppfrag: gammas"<<endl;

 for (kamae_QGS=0; kamae_QGS<=1;kamae_QGS++)
 {
  for (ep= 1.   ;ep<1000.;ep*=epfac)
  for (es= 1.e-2;es<  ep ;es*=esfac)
  {
   reac=0;// pp with Kamae
   result0=ppfrag_spec_int(ep,es,id,reac);
   

  
   reac=1; // pp QGSJET only
   result1    =ppfrag_spec_int(ep,es,id,reac);

 
   // ppfrag v02:  
   result1v02 =ppfrag_v02     (ep,es,id,reac,kamae_QGS); //AWS20130205

   reac=2;// He-p QGSJET only
   result2=ppfrag_spec_int(ep,es,id,reac);

  // ppfrag v02:  interpolated with Kamae
   result2v02 =ppfrag_v02     (ep,es,id,reac,kamae_QGS); //AWS20130205

   reac=3;// p-He QGSJET only
   result3=ppfrag_spec_int(ep,es,id,reac);

  // ppfrag v02:  interpolated with Kamae
   result3v02 =ppfrag_v02     (ep,es,id,reac,kamae_QGS); //AWS20130205

   reac=4;// He-He
  // ppfrag v02:  interpolated with Kamae
   result4v02 =ppfrag_v02     (ep,es,id,reac,kamae_QGS); //AWS20130205

  if(kamae_QGS==0)
  {
   cout <<"test ppfrag v01: ep="<<ep<<" gamma es="<<es
        <<" QGSJET only: pp  =" <<result1
        <<"  p-He="<<result2<<" He-p =" <<result3
        <<endl;
 
   cout <<"test ppfrag v02: model="<<model_name<<" ep="<<ep<<" gamma es="<<es
        <<" QGSJET only: pp  =" <<result1v02
        <<"  p-He="<<result2v02<<" He-p =" <<result3v02<<" He-He=" <<result4v02
	<<"  [GeV, mb/GeV * es(GeV)]"
        <<endl;
   }

  if(kamae_QGS==1)
  {
   cout <<"test ppfrag v01: ep="<<ep<<" gamma es="<<es
        <<" Kamae-QGSJET: pp =" <<result0
        <<endl;
 
   cout <<"test ppfrag v02: model="<<model_name<<" ep="<<ep<<" gamma es="<<es
        <<" Kamae-QGSJET: pp =" <<result1v02
        <<"  p-He="<<result2v02<<" He-p =" <<result3v02<<" He-He=" <<result4v02
	<<"  [GeV, mb/GeV * es(GeV)]"
        <<endl;
   }
 
   cout <<"test ppfrag"<<endl;

   // compare ppfrag with Kamae: should be same below 50 GeV in case of reac=0
	  // Kamae function uses kinetic energy in GeV
          // returns cross-section in mb/GeV times gamma-ray energy in GeV
	 double m_p   = 0.938 ;      // proton mass in GeV  (deprecated hard coding of physics constant !)


        valarray<double> param_a; param_a.resize(9);
        valarray<double> param_b; param_b.resize(8);
        valarray<double> param_c; param_c.resize(5);
        valarray<double> param_d; param_d.resize(5);

         double Tp=ep-m_p; // Kamae notation for kinetic energy of proton in GeV
         param_a=kamae_gamma_param_nd   (Tp);
         param_b=kamae_gamma_param_diff (Tp);
         param_c=kamae_gamma_param_delta(Tp);
         param_d=kamae_gamma_param_res  (Tp);

         int particle=0; // gammas
         double s_nd    = kamae_nd   (es,Tp,param_a,particle);  
         double s_diff  = kamae_diff (es,Tp,param_b);          
         double s_delta = kamae_delta(es,Tp,param_c);          
         double s_res   = kamae_res  (es,Tp,param_d);         
   
                
        // use this notation for consistency with rest of Kachelriess section. 
         double s_kamae= (s_nd + s_diff + s_delta + s_res);

	 // all results here times es

	 // this shows that pp_frag reac=0 uses  (s_nd+s_delta+s_res) for ep < 50 GeV

   cout <<"ppfrag: ep="<<ep<<" es="<<es
        <<"  pp with Kamae="  <<result0
        <<" Kamae routine: Tp="<<Tp
        <<" s_nd="<<s_nd<<" s_diff="<<s_diff<<" s_delta="<<s_delta<<" s_res="<<s_res<<" s_kamae="<<s_kamae
        <<" pp with Kamae / s_nd ="   <<result0/s_nd
        <<" pp with Kamae / (s_nd+s_delta+s_res) ="<<result0/(s_nd+s_delta+s_res)
        <<" pp with Kamae / s_kamae ="<<result0/s_kamae
        <<endl;




  }
 }//kamae_QGSJET

 cout<<"=========== ppfrag: pbar  "<<endl;

  id=1;
  for (ep= 1.;ep<100.;ep*=1.5)
  for (es= 1.;es< ep ;es*=1.5)
  {
   
   reac=0;// pp (no result)     )
   result0=ppfrag_spec_int(ep,es,id,reac);  

    // pbar reac=1,2,3 only
  
   reac=1; // p-p
   result1=ppfrag_spec_int(ep,es,id,reac);

   reac=2;// p-He
   result2=ppfrag_spec_int(ep,es,id,reac);

   reac=3;//He-p
   result3=ppfrag_spec_int(ep,es,id,reac);
   

   cout <<"ppfrag: ep="<<ep<<"  pbar es="<<es
        <<"  pp="  <<result0<<"  p-p =" <<result1
        <<"  p-He="<<result2<<" He-p =" <<result3
       <<endl;
  }


  return 0;
}







//////////////////////////////////////////////////////////////////







//////////////////////////////////////////////////////////////////

int ProductionMatrices::test()
{
  cout<<"ProductionMatrices::test"<<endl;

  string energy_type;
  string energy_units_="MeV";

  energy_type="kinetic_energy";

  init(energy_type,energy_units_);

  energy_type="total_energy";
  init(energy_type,energy_units_);



  /* invalid type: tests check and exit
  energy_type="junk";
  init(energy_type,energy_units_);
  */

  string input_directory="/afs/ipp/u/aws/propagate/c/FITS"; // local test

  input_directory="."; // this is also default if absent in call to read().
  read(input_directory);

  double E_proton_min,E_proton_max, E_proton_factor;
  double E_gamma_min,E_gamma_max, E_gamma_factor;

  E_proton_min   =  1.0;  //  in GeV/nucleon NB total not Ek
  E_proton_min   =  0.4;  //  for case of kinetic energy, but should work also for total
  E_proton_max   =  1.0e3;// = 1   TeV/nucleon
  E_proton_factor=  1.10;

  E_gamma_min    =  1.0e-1;//  in GeV
  E_gamma_max    =  1.0e3; // = 1 TeV
  E_gamma_factor =  1.10;



  int nE_proton_rebin_;
  int nE_gamma_rebin_ ;      

  nE_proton_rebin_ = log(E_proton_max/E_proton_min)/log(E_proton_factor) + 1.001;
  nE_gamma_rebin_  = log(E_gamma_max /E_gamma_min) /log(E_gamma_factor ) + 1.001;

  valarray<double> E_proton_rebin_;
  valarray<double> E_gamma_rebin_ ;

  E_proton_rebin_.resize(nE_proton_rebin_);
  E_gamma_rebin_ .resize( nE_gamma_rebin_);

  for(int i=0;i<nE_proton_rebin_;i++)E_proton_rebin_[i]=  E_proton_min *  pow(E_proton_factor,i);
  for(int i=0;i<nE_gamma_rebin_ ;i++)E_gamma_rebin_ [i]=  E_gamma_min  *  pow(E_gamma_factor ,i);


  int test_rebin =1; // need for next steps !
  if (test_rebin==1)
  { 
   rebin( E_proton_rebin_,  E_gamma_rebin_);
   print(); 
   // return 0; 
  }

  // rebin with energy dispersion
  int test_rebin_with_energy_dispersion  =0;
  if (test_rebin_with_energy_dispersion ==1)
  {
  valarray<double> E_gamma_true;
  int nE_gamma_true=20;
  E_gamma_true.resize( nE_gamma_true);

  for(int i=0;i<nE_gamma_true;i++) E_gamma_true[i]=  E_gamma_min  *  pow(E_gamma_factor ,i);
  rebin_with_energy_dispersion(E_proton_rebin_, E_gamma_rebin_, E_gamma_true  );


  print(); 

  return 0; 
  }


 cout<<" after position G               "<<endl;


  // rebin( E_proton_rebin_,  E_gamma_rebin_,"test units" ); causes exit as required


  valarray<double> proton_spectrum, helium_spectrum;
  proton_spectrum.resize( nE_proton_rebin_);
  helium_spectrum.resize( nE_proton_rebin_);

  int             nE_gamma_band;
  valarray<double> E_gamma_band_min_;
  valarray<double> E_gamma_band_max_;


  valarray<double> emissivity_total_;
  emissivity_total_.resize( nE_gamma_rebin_);
 

  // Shikaze et al 2007 http://adsabs.harvard.edu/abs/2007APh....28..154S
  // Fig 8 and 9 parameters
  // flux = A beta^P1 R^-P2

  cout<<"emissivity from Shikaze etal protons + Helium"<<endl;

  // protons
  double A=1.94e4*1e-4;// from m-2 sr-1 s-1 GeV-1 to cm-2  sr-1 s-1 GeV-1
  double P1=0.70;
  double P2=2.76;
  double m_p=0.938;// proton mass in GeV


 int test_without_CR_Spectrum=0;
 if( test_without_CR_Spectrum==1)
 {
  proton_spectrum=0.;
  helium_spectrum=0.;

  for(int i=0;i<nE_proton_rebin_;i++)
    {
      double Ek   = E_proton_rebin_[i]-m_p;// kinetic energy per nucleon
      double g    = E_proton_rebin_[i]/m_p;// Lorentz factor of proton
      double beta = sqrt(1.0-1.0/g/g);                  // velocity/c     of proton
      double R    = sqrt(pow(E_proton_rebin_[i],2)-pow(m_p,2)); // rigidity = momentum

      proton_spectrum[i]=A*pow(beta,P1)*pow(R,-P2);

      

     cout<<"E total= "<< E_proton_rebin_[i]<< " GeV Ek="<<Ek<<" GeV beta="<<beta<<" R="<<R 
                      <<" Shikaze proton spectrum="<< proton_spectrum[i]<<"  cm-2  sr-1 s-1 GeV-1"<<endl;

    }


  // Helium

         A=7.10e3*1e-4;// from m-2  sr-1 s-1 (GeV/n)-1 to cm-2  sr-1 s-1 (GeV/n)-1
         P1=0.50;
         P2=2.78;

  for(int i=0;i<nE_proton_rebin_;i++)
    {
      double Ek, g,beta,R;

      
	{

	     double Z=2.0;
             Ek   = E_proton_rebin_[i] - m_p;                    // kinetic energy per nucleon
	     g    = E_proton_rebin_[i]/m_p;                   // Lorentz factor of helium
             beta = sqrt(1.0-1.0/g/g);                               // velocity/c     of helium
             R    = sqrt(pow(4.*E_proton_rebin_[i],2)-pow(4.*m_p,2))/Z; // rigidity = particle momentum/Z  AWS20120319

             helium_spectrum[i]=A*pow(beta,P1)  *pow(R,-P2) ; 
	}
      cout<<"E total= "<< E_proton_rebin_[i]<< " GeV Ek/n="<<Ek  <<" GeV beta="<<beta<<" R="<<R 
                       <<" Shikaze helium spectrum="<< helium_spectrum[i]<<"  cm-2  sr-1 s-1 (GeV/n)-1"<<endl;
    }

  

 

  
  
  //  for(int i=0;i<nE_proton_rebin_;i++)cout<<"actually used for emissivity: E total= "<< E_proton_rebin_[i]
  //                                           <<" protons="<<proton_spectrum[i]<<" helium="<<helium_spectrum[i]<<endl;


    emissivity                         (proton_spectrum,helium_spectrum); // p and He separately

    emissivity_total_ =emissivity_total(proton_spectrum,helium_spectrum); // sum of p and He

    cout<<endl<<"after first rebin and emissivity"<<endl;
    print();

  }// test_without_CR_Spectrum

 cout<<" after position F               "<<endl;

    // define energy bands

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
  


  // Table 1 of Abdo etal 703, 1249, 2009 (HI emissivity paper)
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
   emiss_measured[i]*=1.0e-3; // MeV-> GeV
   emiss_measured[i]*=(E_gamma_band_max_[i]- E_gamma_band_min_[i]);
   emiss_measured[i]/=(E_gamma_band_max_[i]* E_gamma_band_min_[i]);
   cout<<"Abdo etal 2009 703, 1249, 2009 (HI emissivity paper) measured emissivity: E_gamma="
        <<E_gamma_band_min_[i]<<" - "<<E_gamma_band_max_[i]
        <<" GeV =" << emiss_measured[i]<<endl;
  }

  // Table 2 of Ackermann etal ApJ 726, 81, 2011   (3rd quad paper)
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
   emiss_measured[i]*=1.0e-3; // MeV-> GeV
   emiss_measured[i]*=(E_gamma_band_max_[i]- E_gamma_band_min_[i]);
   emiss_measured[i]/=(E_gamma_band_max_[i]* E_gamma_band_min_[i]);
   cout<<"Ackermann etal ApJ 726, 81, 2011   (3rd quad paper) measured emissivity: E_gamma="
        <<E_gamma_band_min_[i]<<" - "<<E_gamma_band_max_[i]
        <<" GeV =" << emiss_measured[i]<<endl;
  }

  // finished defining energies and bands

 cout<<" after position E               "<<endl;


  double E_gamma_band_min;
  double E_gamma_band_max;
  double E_gamma_band_factor=1.01;// for integration over gamma-ray band

  double emiss_band, emiss_band_protons, emiss_band_helium;

  valarray<double>zero_spectrum;
  zero_spectrum.resize(proton_spectrum.size());
  zero_spectrum=0.;


  int  test_proton_helium_separately=0;

  if ( test_proton_helium_separately==1)
  {



  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
    {
     E_gamma_band_min=E_gamma_band_min_[iE_gamma_band];
     E_gamma_band_max=E_gamma_band_max_[iE_gamma_band];

         emiss_band=
    emissivity_band( proton_spectrum,  helium_spectrum,
                          E_proton_rebin_,
                          E_gamma_band_min, E_gamma_band_max, E_gamma_band_factor);

        emiss_band_protons=
    emissivity_band( proton_spectrum,  zero_spectrum,
                          E_proton_rebin_,
                          E_gamma_band_min, E_gamma_band_max, E_gamma_band_factor);

        emiss_band_helium=
    emissivity_band( zero_spectrum,  helium_spectrum,
                          E_proton_rebin_,
                          E_gamma_band_min, E_gamma_band_max, E_gamma_band_factor);


	 cout<<"test: Shikaze LIS "                 
             <<" E_gamma_band_min="<< E_gamma_band_min<<" E_gamma_band_max="<< E_gamma_band_max
             <<" emiss_band="<<emiss_band      <<" emiss_band_protons="<<emiss_band_protons    
             <<" emiss_band_helium="<<emiss_band_helium         
             <<" predicted/measured = " <<emiss_band/emiss_measured[iE_gamma_band]
             <<endl;

	 cout<<"Shikaze LIS "                 
             <<" E_gamma="<< E_gamma_band_min<<" -"<< E_gamma_band_max
             <<" GeV predicted emiss="  <<emiss_band   
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band/emiss_measured[iE_gamma_band]
             <<endl;

}


  }//test_proton_helium_separately



  // test multiple bands version

  valarray<double> emiss_band_;
  emiss_band_.resize(nE_gamma_band);

  int  test_multiple_bands=0;

  if(test_multiple_bands==1)
  {

  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum,
                          E_proton_rebin_,
                          E_gamma_band_min_, E_gamma_band_max_, E_gamma_band_factor);
  

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" output from multiple bands version: Shikaze LIS "                 
             <<" E_gamma="<< E_gamma_band_min_[iE_gamma_band]<<" -"<< E_gamma_band_max_[iE_gamma_band]
             <<" GeV predicted emiss="  <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;

  }//test_multiple_bands==1

 cout<<" after position D               "<<endl;


 ////////////////////////////////////////////////


  int test_varying_modulation=0;

 if (test_varying_modulation==1)
 {

 cout<<"emissivity from varying modulation of  protons + Helium"<<endl;

 /*
Gleeson & Axford,  ApJ 154, 1011 (1968)
The formulae here are based on equation (16) of that paper.

 J(E)           J(E + |Z| e phi)
________   =    _________________________ 
E^2-Eo^2         (E + |Z| e phi)^2 - Eo^2        

where E = total energy of particle, Eo = rest mass energy of particle.

see also galplot: modulate.cc
 */

 // Shikaze LIS uses modulation 591 MV for 1998 data, so compute relative to this.

 for(double Phi=0.;Phi< 600.;Phi+= 50.)
 {

   //          Phi= 91. ;// modulation parameter in standard MV units
  double PhiGV=Phi/1000.;// modulation parameter in GV 

   A=1.94e4*1e-4;// from m-2 sr-1 s-1 GeV-1 to cm-2  sr-1 s-1 GeV-1
   P1=0.70;
   P2=2.76;
   m_p=0.938;// proton mass in GeV

   proton_spectrum=0.;
   helium_spectrum=0.;

  for(int i=0;i<nE_proton_rebin_;i++)
    {
      double Ek   = E_proton_rebin_[i]-m_p;// kinetic energy per nucleon
      double g    = E_proton_rebin_[i]/m_p;// Lorentz factor of proton
      double beta = sqrt(1.0-1.0/g/g);                  // velocity/c     of proton
      double R    = sqrt(pow(E_proton_rebin_[i],2)-pow(m_p,2)); // rigidity = momentum

      double        E_prime = E_proton_rebin_[i]+PhiGV;// total energy prime
      double        g_prime = E_prime/m_p;
      double     beta_prime = sqrt(1.0-1.0/g_prime/g_prime);   
      double        R_prime = sqrt(pow(E_prime,2)-pow(m_p,2)); // rigidity = momentum
      double        factor  = ( pow(  E_proton_rebin_[i],2)-pow(m_p,2)) / ( pow(  E_prime,2)-pow(m_p,2))   ;

      proton_spectrum[i]=A*pow(beta_prime,P1)*pow(R_prime,-P2)* factor;

      cout<<"E total = "<< E_proton_rebin_[i] <<" Ek="<<Ek<< " E_prime="<<E_prime<<" R_prime="<<R_prime
                        <<" Shikaze modulated proton spectrum="<< proton_spectrum[i]<<"  cm-2  sr-1 s-1 GeV-1" <<endl;
    }

  // Helium

         A=7.10e3*1e-4;// from m-2  sr-1 s-1 GeV-1 to cm-2  sr-1 s-1 GeV-1
         P1=0.50;
         P2=2.78;

  for(int i=0;i<nE_proton_rebin_;i++)
    {
     
     
      {
      double Z=2.;
      double Ek   = E_proton_rebin_[i] - m_p;// kinetic energy per nucleon
      double g    = E_proton_rebin_[i] / m_p;// Lorentz factor of particle
      double beta = sqrt(1.0-1.0/g/g);                  // velocity/c     of particle
      double R    = sqrt(pow(4.*E_proton_rebin_[i],2)-pow(4.*m_p,2))/Z;// rigidity = particle momentum/Z 
    


      double        E_prime = E_proton_rebin_[i] +Z *PhiGV/4.0;        // total energy per nucleon prime
      double        g_prime = E_prime/m_p;                             // Lorentz factor of particle
      double     beta_prime = sqrt(1.0-1.0/g_prime/g_prime);           // velocity/c     of particle  
      double        R_prime = sqrt(pow(4.*E_prime,2)-pow(4.*m_p,2))/Z; // rigidity = particle momentum/Z
      double        factor  = ( pow(  4.*E_proton_rebin_[i],2)-pow(4.*m_p,2)) / ( pow(  4.*E_prime,2)-pow(4.*m_p,2))   ;

      helium_spectrum[i]=A*pow(beta_prime,P1)*pow(R_prime,-P2)* factor;

      cout<<"E total = "<< E_proton_rebin_[i] <<" Ek="<<Ek<<" GeV/n  E_prime="<<E_prime<<" R_prime="<<R_prime
                        <<" Shikaze modulated helium spectrum="<< helium_spectrum[i]<<"  cm-2  sr-1 s-1 GeV-1" <<endl;
      }


     
    }

  
  
  

         

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
    {
     E_gamma_band_min=E_gamma_band_min_[iE_gamma_band];
     E_gamma_band_max=E_gamma_band_max_[iE_gamma_band];

         emiss_band=
    emissivity_band( proton_spectrum,  helium_spectrum,
                          E_proton_rebin_,
                          E_gamma_band_min, E_gamma_band_max, E_gamma_band_factor);

        emiss_band_protons=
    emissivity_band( proton_spectrum,  zero_spectrum,
                          E_proton_rebin_,
                          E_gamma_band_min, E_gamma_band_max, E_gamma_band_factor);

        emiss_band_helium=
    emissivity_band( zero_spectrum,  helium_spectrum,
                          E_proton_rebin_,
                          E_gamma_band_min, E_gamma_band_max, E_gamma_band_factor);


	cout<<"test: Shikaze modulated Phi="<<Phi<<" delta Phi="<<591-Phi
             <<" E_gamma_band_min="<< E_gamma_band_min<<" E_gamma_band_max="<< E_gamma_band_max
             <<" emiss_band="<<emiss_band      <<" emiss_band_protons="<<emiss_band_protons    
             <<" emiss_band_helium="<<emiss_band_helium                        <<endl;



	 cout<<"Shikaze modulated Phi="<<Phi<<" delta Phi="<<591-Phi
             <<" E_gamma="<< E_gamma_band_min<<" -"<< E_gamma_band_max
             <<" GeV predicted emiss="  <<emiss_band   
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band/emiss_measured[iE_gamma_band]
             <<endl;
}

  // test multiple band version

  emiss_band_.resize(nE_gamma_band);
  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum,
                          E_proton_rebin_,
                          E_gamma_band_min_, E_gamma_band_max_, E_gamma_band_factor);
  

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<"result from multiple bands version: Shikaze modulated Phi="<<Phi<<" delta Phi="<<591-Phi         
             <<" E_gamma="<< E_gamma_band_min_[iE_gamma_band]<<" -"<< E_gamma_band_max_[iE_gamma_band]
             <<" GeV predicted emiss="  <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;
  

 }// for Phi


 }//test_varying_modulation==1

  ////////////////////////////////////////////////


 cout<<" after position C               "<<endl;



 // test CR_Spectrum

  CR_Spectrum cr_spectrum;
  cr_spectrum.test();

  cout<< "--- ProductionMatrices:: test : using CR_Spectrum"<<endl;

  

  string name;
 
  valarray<double>parameters;
  parameters.resize(10);

  //////////////////////////////////////////////////////

  int test_using_CR_spectrum=0;

  if(test_using_CR_spectrum==1)
  {

 for(double Phi=300.;Phi< 620.;Phi+= 20.)
 {

  parameters[0]=Phi ;// modulation in MV

  int debug=1;

  name="Shikaze protons";
  cr_spectrum.init(name,parameters, debug);
  proton_spectrum=cr_spectrum.flux(E_proton_rebin);

  name="Shikaze Helium";
  cr_spectrum.init(name,parameters, debug);
  helium_spectrum=cr_spectrum.flux(E_proton_rebin);

  emiss_band_.resize(nE_gamma_band);
  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum,
                          E_proton_rebin_,
                          E_gamma_band_min_, E_gamma_band_max_, E_gamma_band_factor);
  

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" result from multiple bands CR_Spectrum  Shikaze modulated"   
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_[iE_gamma_band]<<" -"<< E_gamma_band_max_[iE_gamma_band]
             <<" GeV predicted emiss="  <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;

 

  }// Phi


  }//test_using_CR_spectrum==1

 cout<<" after position B               "<<endl;

  /////////////////////////////////////

// test MeV version

  cout<<"testing MeV version"<<endl;

  energy_units_="MeV";
  valarray<double>  E_proton_rebin_MeV;
  E_proton_rebin_MeV.resize(E_proton_rebin.size());
  E_proton_rebin_MeV      = E_proton_rebin*1.e3; //GeV->MeV

  valarray<double>  E_gamma_rebin_MeV;
  E_gamma_rebin_MeV.resize(E_gamma_rebin.size());
  E_gamma_rebin_MeV      = E_gamma_rebin*1.e3; //GeV->MeV


  valarray<double> E_gamma_band_min_MeV,E_gamma_band_max_MeV;
  E_gamma_band_min_MeV.resize(nE_gamma_band);
  E_gamma_band_max_MeV.resize(nE_gamma_band);
  E_gamma_band_min_MeV =E_gamma_band_min_*1e3; //GeV->MeV
  E_gamma_band_max_MeV =E_gamma_band_max_*1e3; //GeV->MeV


  double Phi=591;// test case of no modulation difference first
  parameters[0]=Phi ;// modulation in MV



  name="Shikaze protons";
  int debug=0;
  cr_spectrum.init(name,parameters, debug);
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_); // still GeV but can use like this. use E_proton_rebin_ not E_proton_rebin since latter is in class
  
  name="Shikaze Helium";
  cr_spectrum.init(name,parameters, debug);
  helium_spectrum=cr_spectrum.flux(E_proton_rebin_); // still GeV but can use like this

  // GeV units
  energy_units_="GeV";

  /* test rebin and emissivity
  rebin( E_proton_rebin , E_gamma_rebin);

  emissivity(proton_spectrum, helium_spectrum);
  cout<<"checking matrices after rebin and emissivity for GeV"<<endl;
  print();
  */


 for(int i=0;i<nE_proton_rebin_;i++)  cout<<"actually used for checking  emissivity: E_proton_rebin "
				 	  << E_proton_rebin[i]<<" E_proton_rebin_=  "<< E_proton_rebin_[i]
                                          <<" E_proton_rebin_MeV=  "<< E_proton_rebin_MeV[i]
                                          <<" protons="<<proton_spectrum[i]<<" helium="<<helium_spectrum[i]
                                          <<" energy_units="<<energy_units_  <<endl;

  emiss_band_.resize(nE_gamma_band);
  emiss_band_=0.;
  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum,
                                  E_proton_rebin,
				  E_gamma_band_min_, E_gamma_band_max_, E_gamma_band_factor);
  
  

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" checking result for "  
	     <<" energy_units="<<energy_units 
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_[iE_gamma_band]<<" -"<< E_gamma_band_max_[iE_gamma_band]
             <<" predicted emiss="  <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;
 
 cout<<" after position A               "<<endl;


  /////////////////////////////////////////////

  // MeV units

  /////////////////////////////////////////////

  energy_units_="MeV";
  proton_spectrum*=1.e-3;//GeV^-1 -> MeV^-1
  helium_spectrum*=1.e-3;//GeV^-1 -> MeV^-1

  /* test rebin and emissivity
  rebin( E_proton_rebin_MeV , E_gamma_rebin_MeV);

  emissivity(proton_spectrum, helium_spectrum);
  cout<<"checking after matrices after rebin and emissivity for MeV"<<endl;

  print();
  */


  for(int i=0;i<nE_proton_rebin_;i++)cout<<"actually used for checking  emissivity: E_proton_rebin "
					 << E_proton_rebin[i]<<" E_proton_rebin_=  "<< E_proton_rebin_[i]
                                          <<" E_proton_rebin_MeV=  "<< E_proton_rebin_MeV[i]
                                          <<" protons="<<proton_spectrum[i]<<" helium="<<helium_spectrum[i]
                                          <<" energy_units="<<energy_units_  <<endl;
  emiss_band_.resize(nE_gamma_band);
  emiss_band_=0.;
  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum,
                                  E_proton_rebin_MeV,
				   E_gamma_band_min_MeV, E_gamma_band_max_MeV, E_gamma_band_factor);
  

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" checking result for "  
	     <<" energy_units="<<energy_units 
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="  <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;

  // the above comparison showed correct emissivities also in bands.

  // for(double Phi=300.;Phi< 620.;Phi+= 20.)
   Phi=591;// test case of no modulation difference first
 {

  parameters[0]=Phi ;// modulation in MV

  debug=1;

  name="Shikaze protons";
  cr_spectrum.init(name,parameters, debug);
  //  proton_spectrum=cr_spectrum.flux(E_proton_rebin); // still GeV but can use after conversion
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_); // NB previous invocation of emissivity_band has modified E_proton_bin in class
  proton_spectrum*=1.e-3;//GeV^-1 -> MeV^-1

  name="Shikaze Helium";
  cr_spectrum.init(name,parameters, debug);
  //  helium_spectrum=cr_spectrum.flux(E_proton_rebin); // still GeV but can use after conversion
  helium_spectrum=cr_spectrum.flux(E_proton_rebin_); // NB previous invocation of emissivity_band has modified E_proton_bin in class
  helium_spectrum*=1.e-3;//GeV^-1 -> MeV^-1

  energy_units_="MeV";

  

  for(int i=0;i<nE_proton_rebin_;i++)cout<<"actually used for checking again emissivity: E_proton_rebin "
					 << E_proton_rebin[i]<<" E_proton_rebin_=  "<< E_proton_rebin_[i]
                                          <<" E_proton_rebin_MeV=  "<< E_proton_rebin_MeV[i]
                                          <<" protons="<<proton_spectrum[i]<<" helium="<<helium_spectrum[i]
                                          <<" energy_units="<<energy_units_  <<endl;


  emiss_band_.resize(nE_gamma_band);
  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum,
                                  E_proton_rebin_MeV,
				   E_gamma_band_min_MeV, E_gamma_band_max_MeV, E_gamma_band_factor);
  

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" checking again result from multiple for "  
	     <<" energy_units="<<energy_units_ 
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="  <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;

 

  }// Phi


 // Shikaze CR_Spectrum MeV version

  Phi=591;// test case of no modulation difference first

 for(double Phi=300.;Phi< 620.;Phi+= 20.)
 {

  parameters[0]=Phi ;// modulation in MV

  debug=1;

  name="Shikaze protons MeV";
  cr_spectrum.init(name,parameters, debug);
  
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 
  

  name="Shikaze Helium MeV";
  cr_spectrum.init(name,parameters, debug);

  helium_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 
  

  energy_units_="MeV";

  
  for(int i=0;i<nE_proton_rebin_;i++)cout<<"actually used for checking Shikaze MeV CR_Spectrum emissivity: E_proton_rebin= "
					 << E_proton_rebin[i]<<" E_proton_rebin_=  "<< E_proton_rebin_[i]
                                          <<" E_proton_rebin_MeV=  "<< E_proton_rebin_MeV[i]
					  <<" Phi="<<Phi
                                          <<" protons="<<proton_spectrum[i]<<" helium="<<helium_spectrum[i]
                                          <<" energy_units="<<energy_units_  <<endl;


  emiss_band_.resize(nE_gamma_band);
  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum,
                                  E_proton_rebin_MeV,
				   E_gamma_band_min_MeV, E_gamma_band_max_MeV, E_gamma_band_factor);
  

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" checking Shikaze MeV CR_Spectrum for "  
	     <<" energy_units="<<energy_units_ 
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="  <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;

 

  }// Phi


 cout<<"after  Phi                      "<<endl;

 ///////////////////////////////////////////

 int test_gamma_bands_matrices =0;

 if (test_gamma_bands_matrices==1)
 {
   cout<<"testing emissivity for  gamma-ray energy bands using precomputed matrices"<<endl;

   cout<<"case1"<<endl;

   /* to avoid setting directly
  energy_units="GeV";

  gen_gamma_band_matrices(E_proton_rebin_,
                          E_gamma_band_min_, E_gamma_band_max_,  E_gamma_band_factor);
   */

   // this should not work since init(energy_type,energy_units) has not been called: test this 
  
  energy_units_="GeV";

  
  gen_gamma_band_matrices(E_proton_rebin_,
                          E_gamma_band_min_, E_gamma_band_max_,  E_gamma_band_factor);
  
   

  Phi=591;
  parameters[0]=Phi ;// modulation in MV



  name="Shikaze protons";
  int debug=0;
  cr_spectrum.init(name,parameters, debug);
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_); // GeV
  
  name="Shikaze Helium";
  cr_spectrum.init(name,parameters, debug);
  helium_spectrum=cr_spectrum.flux(E_proton_rebin_); 

  emiss_band_.resize(nE_gamma_band);
  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum);


  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" using precomputed matrices "  
	     <<" energy_units="<<energy_units 
	     <<" energy_units_="<<energy_units_ 
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_[iE_gamma_band]<<" -"<< E_gamma_band_max_[iE_gamma_band]
             <<" predicted emiss="      <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;

  /////////////////////

    cout<<"case2"<<endl;

  energy_units="MeV";

  gen_gamma_band_matrices(E_proton_rebin_MeV,
                          E_gamma_band_min_MeV, E_gamma_band_max_MeV,  E_gamma_band_factor);

  name="Shikaze protons MeV";
  cr_spectrum.init(name,parameters, debug);
  
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 
  

  name="Shikaze Helium MeV";
  cr_spectrum.init(name,parameters, debug);

  helium_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 

  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum);

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" using precomputed matrices "  
	     <<" energy_units="<<energy_units 
	     <<" energy_type ="<<energy_type
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="      <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;


  ///////////////
    cout<<"case3"<<endl;

  energy_units="MeV";

  energy_type="kinetic_energy";
  init(energy_type,energy_units);

  

  gen_gamma_band_matrices(E_proton_rebin_MeV,
                          E_gamma_band_min_MeV, E_gamma_band_max_MeV,  E_gamma_band_factor);

  name="Shikaze protons kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);
  
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 
  

  name="Shikaze Helium kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);

  helium_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 

  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum);

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" using precomputed matrices " 
 	     <<" model_name="<<model_name
	     <<" energy_units="<<energy_units 
	     <<" energy_type ="<<energy_type  
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="      <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;


 /////////////// Kachelriess model
  int Kachelriess =0;
 if ( Kachelriess==1)
 {
  energy_units="MeV";

  energy_type="kinetic_energy";
  string model_name="Kachelriess";
  init(energy_type,energy_units,0,model_name);

  

  gen_gamma_band_matrices(E_proton_rebin_MeV,
                          E_gamma_band_min_MeV, E_gamma_band_max_MeV,  E_gamma_band_factor);

  name="Shikaze protons kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);
  
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 
  

  name="Shikaze Helium kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);

  helium_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 

  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum);

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" using precomputed matrices " 
	     <<" model_name="<<model_name
	     <<" energy_units="<<energy_units 
	     <<" energy_type ="<<energy_type  
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="      <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;
 }
 /////////////// Kachelriess2 model

 int Kachelriess2=0;
 if(Kachelriess2==1)
 {
  energy_units="MeV";

  energy_type="kinetic_energy";
  model_name="Kachelriess2";
  init(energy_type,energy_units,0,model_name);

  

  gen_gamma_band_matrices(E_proton_rebin_MeV,
                          E_gamma_band_min_MeV, E_gamma_band_max_MeV,  E_gamma_band_factor);

  name="Shikaze protons kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);
  
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 
  

  name="Shikaze Helium kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);

  helium_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 

  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum);

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" using precomputed matrices " 
	     <<" model_name="<<model_name
	     <<" energy_units="<<energy_units 
	     <<" energy_type ="<<energy_type  
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="      <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;


 }

 /////////////// Dermer-Kachelriess model
    int Dermer_Kachelriess=0;

  if( Dermer_Kachelriess==1)
  {
  energy_units="MeV";

  energy_type="kinetic_energy";
  model_name="Dermer_Kachelriess";
  init(energy_type,energy_units,0,model_name);

  

  gen_gamma_band_matrices(E_proton_rebin_MeV,
                          E_gamma_band_min_MeV, E_gamma_band_max_MeV,  E_gamma_band_factor);

  name="Shikaze protons kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);
  
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 
  

  name="Shikaze Helium kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);

  helium_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 

  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum);

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" using precomputed matrices " 
	     <<" model_name="<<model_name
	     <<" energy_units="<<energy_units 
	     <<" energy_type ="<<energy_type  
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="      <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;

 }// if  Dermer_Kachelriess==1


/////////////// Kachelriess2013 model
    int Kachelriess2013=1;

    if(Kachelriess2013==1)
  {
    cout<<"Kachelriess2013"<<endl;

  energy_units="MeV";

  energy_type="kinetic_energy";
  model_name="Kachelriess2013";
  init(energy_type,energy_units,0,model_name);

  

  gen_gamma_band_matrices(E_proton_rebin_MeV,
                          E_gamma_band_min_MeV, E_gamma_band_max_MeV,  E_gamma_band_factor);

  name="Shikaze protons kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);
  
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 
  

  name="Shikaze Helium kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);

  helium_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 

  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum);

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" using precomputed matrices " 
	     <<" model_name="<<model_name
	     <<" energy_units="<<energy_units 
	     <<" energy_type ="<<energy_type  
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="      <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;

 }// if  Kachelriess2013==1










 /////////////// FLUKA model
 int FLUKA=0;
 if(FLUKA==1)
 {
   cout<<"testing FLUKA"<<endl;
  energy_units="MeV";

  energy_type="kinetic_energy";
  model_name="FLUKA";
  init(energy_type,energy_units,0,model_name);

  

  gen_gamma_band_matrices(E_proton_rebin_MeV,
                          E_gamma_band_min_MeV, E_gamma_band_max_MeV,  E_gamma_band_factor);

  name="Shikaze protons kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);
  
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 
  

  name="Shikaze Helium kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);

  helium_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 

  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum);

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" using precomputed matrices " 
	     <<" model_name="<<model_name
	     <<" energy_units="<<energy_units 
	     <<" energy_type ="<<energy_type  
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="      <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;

 }//FLUKA

 }//test_gamma_bands_matrices

 cout<<"after  test_gamma_bands_matrices"<<endl;

 ///////////////////////////////////////////////////////
   int test_Gaisser_protons=0;

 if (test_Gaisser_protons==1)
 {

 for(double Phi=0;Phi<=500.;Phi+= 50.)
 {

  parameters[0]=Phi ;// modulation in MV

  int debug=1;

  name="Gaisser protons";
  cr_spectrum.init(name,parameters, debug);
  proton_spectrum=cr_spectrum.flux(E_proton_rebin);

  name="Gaisser Helium high";
  cr_spectrum.init(name,parameters, debug);
  helium_spectrum=cr_spectrum.flux(E_proton_rebin);


  emiss_band_.resize(nE_gamma_band);
  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum,
                          E_proton_rebin_,
                          E_gamma_band_min_, E_gamma_band_max_, E_gamma_band_factor);
  

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<"CR_Spectrum  Gaisser demodulated, p+He high"   
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_[iE_gamma_band]<<" -"<< E_gamma_band_max_[iE_gamma_band]
             <<" GeV predicted emiss="  <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;

  name="Gaisser Helium low";
  cr_spectrum.init(name,parameters, debug);
  helium_spectrum=cr_spectrum.flux(E_proton_rebin);


  emiss_band_.resize(nE_gamma_band);
  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum,
                          E_proton_rebin_,
                          E_gamma_band_min_, E_gamma_band_max_, E_gamma_band_factor);
  

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<"CR_Spectrum  Gaisser demodulated, p+He low "   
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_[iE_gamma_band]<<" -"<< E_gamma_band_max_[iE_gamma_band]
             <<" GeV predicted emiss="  <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;



  }// Phi


 }// test_Gaisser_protons==1

 ///////////////////////////////////////////////////////////////

  int test_ppfrag_ =0;
  if( test_ppfrag_==1)
  {

 // Kachelreiss model (ppfrag)

  energy_units="MeV";
  energy_type="kinetic_energy";
  string model_name="Kachelriess";
  init(energy_type,energy_units,0,model_name);

  test_ppfrag();


  
  energy_units="MeV";
  energy_type="kinetic_energy";
  model_name="Kachelriess2013";
  init(energy_type,energy_units,0,model_name);

  test_ppfrag();

  }
  /////////////////////////////////////

  // rebin with energy dispersion
      test_rebin_with_energy_dispersion  =0;
  if (test_rebin_with_energy_dispersion ==1)
  {

    cout<<"test rebin_with_energy_dispersion.."<<endl;

  energy_units="MeV";
  energy_type="kinetic_energy";
  model_name="Kachelriess2013";
  init(energy_type,energy_units,0,model_name);


  E_proton_min   =  100;  //  in MeV/nucleon NB total not Ek
  E_proton_max   =  1.0e5;
  E_proton_factor=  1.20;

  E_gamma_min    =  100.;   //  in MeV
  E_gamma_max    =  1.0e5;                 
  E_gamma_factor =  1.05;



 

  nE_proton_rebin_ = log(E_proton_max/E_proton_min)/log(E_proton_factor) + 1.001;
  nE_gamma_rebin_  = log(E_gamma_max /E_gamma_min) /log(E_gamma_factor ) + 1.001;


  E_proton_rebin_.resize(nE_proton_rebin_);
  E_gamma_rebin_ .resize( nE_gamma_rebin_);

  for(int i=0;i<nE_proton_rebin_;i++)E_proton_rebin_[i]=  E_proton_min *  pow(E_proton_factor,i);
  for(int i=0;i<nE_gamma_rebin_ ;i++)E_gamma_rebin_ [i]=  E_gamma_min  *  pow(E_gamma_factor ,i);



  valarray<double> E_gamma_true;
  int nE_gamma_true=150;
  E_gamma_true.resize( nE_gamma_true);
                
              
  E_gamma_factor =  1.05;
  for(int i=0;i<nE_gamma_true;i++) E_gamma_true[i]=  1.0* E_gamma_min  *  pow(E_gamma_factor ,i);
  rebin_with_energy_dispersion(E_proton_rebin_, E_gamma_rebin_, E_gamma_true  );


  //print(); // done in routine

  
  }

// test emissivity bands for energy dispersion
/////////////// Kachelriess2013 model
    int Kachelriess2013=1;

    if(Kachelriess2013==1)
  {
    cout<<"test Kachelriess2013 with energy dispersion"<<endl;

  energy_units="MeV";

  energy_type="kinetic_energy";
  model_name="Kachelriess2013";
  init(energy_type,energy_units,0,model_name);


  E_gamma_band_min_MeV[0]=10;  E_gamma_band_max_MeV[0]=20;
  E_gamma_band_min_MeV[1]=20;  E_gamma_band_max_MeV[1]=30;
  E_gamma_band_min_MeV[2]=30;  E_gamma_band_max_MeV[2]=40;
  E_gamma_band_min_MeV[3]=40;  E_gamma_band_max_MeV[3]=50;
  E_gamma_band_min_MeV[4]=50;  E_gamma_band_max_MeV[4]=60;

  gen_gamma_band_matrices(E_proton_rebin_MeV,
                          E_gamma_band_min_MeV, E_gamma_band_max_MeV,  E_gamma_band_factor);

  name="Shikaze protons kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);
  
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 
  

  name="Shikaze Helium kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);

  helium_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 

  debug=0;// class member

  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum);

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" using precomputed matrices " 
	     <<" model_name="<<model_name
	     <<" energy_dispersion="<<energy_dispersion
	     <<" energy_units="<<energy_units 
	     <<" energy_type ="<<energy_type  
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="      <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;


  valarray<double> emiss_band_nodisp;
  emiss_band_nodisp.resize(emiss_band_.size());
  emiss_band_nodisp=emiss_band_;



  ///

  int energy_dispersion_=1; // above defaults to 1
  init(energy_type,energy_units,0,model_name,energy_dispersion_);



  gen_gamma_band_matrices(E_proton_rebin_MeV,
                          E_gamma_band_min_MeV, E_gamma_band_max_MeV,  E_gamma_band_factor);

  name="Shikaze protons kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);
  
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 
  

  name="Shikaze Helium kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);

  helium_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 

  debug=0;// class member

  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum);

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" using precomputed matrices " 
	     <<" model_name="<<model_name
	     <<" energy_dispersion="<<energy_dispersion
	     <<" energy_units="<<energy_units 
	     <<" energy_type ="<<energy_type  
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="      <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
	     <<" edisp/noedisp ="       <<emiss_band_[iE_gamma_band] /emiss_band_nodisp[iE_gamma_band] 
             <<endl;















 }// if  Kachelriess2013==1









/////////////// Dermer_Kachelriess2013 model
   int Dermer_Kachelriess2013 =0;

    if(Dermer_Kachelriess2013==1)
  {
    cout<<"test Dermer_Kachelriess2013 with energy dispersion"<<endl;

  energy_units="MeV";

  energy_type="kinetic_energy";
  model_name="Dermer_Kachelriess2013";
  init(energy_type,energy_units,0,model_name);


  E_gamma_band_min_MeV[0]=10;  E_gamma_band_max_MeV[0]=20;
  E_gamma_band_min_MeV[1]=20;  E_gamma_band_max_MeV[1]=30;
  E_gamma_band_min_MeV[2]=30;  E_gamma_band_max_MeV[2]=40;
  E_gamma_band_min_MeV[3]=40;  E_gamma_band_max_MeV[3]=50;
  E_gamma_band_min_MeV[4]=50;  E_gamma_band_max_MeV[4]=60;

  gen_gamma_band_matrices(E_proton_rebin_MeV,
                          E_gamma_band_min_MeV, E_gamma_band_max_MeV,  E_gamma_band_factor);

  name="Shikaze protons kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);
  
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 
  

  name="Shikaze Helium kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);

  helium_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 

  debug=0;// class member

  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum);

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" using precomputed matrices " 
	     <<" model_name="<<model_name
	     <<" energy_dispersion="<<energy_dispersion
	     <<" energy_units="<<energy_units 
	     <<" energy_type ="<<energy_type  
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="      <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;


  valarray<double> emiss_band_nodisp;
  emiss_band_nodisp.resize(emiss_band_.size());
  emiss_band_nodisp=emiss_band_;



  ///

  int energy_dispersion_=1; // above defaults to 1
  init(energy_type,energy_units,0,model_name,energy_dispersion_);



  gen_gamma_band_matrices(E_proton_rebin_MeV,
                          E_gamma_band_min_MeV, E_gamma_band_max_MeV,  E_gamma_band_factor);

  name="Shikaze protons kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);
  
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 
  

  name="Shikaze Helium kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);

  helium_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 

  debug=0;// class member

  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum);

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" using precomputed matrices " 
	     <<" model_name="<<model_name
	     <<" energy_dispersion="<<energy_dispersion
	     <<" energy_units="<<energy_units 
	     <<" energy_type ="<<energy_type  
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="      <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
	     <<" edisp/noedisp ="       <<emiss_band_[iE_gamma_band] /emiss_band_nodisp[iE_gamma_band] 
             <<endl;















 }// if  Dermer_Kachelriess2013==1

/////////////// Dermer_Kachelriess2013_heavy model
   int Dermer_Kachelriess2013_heavy =1;

    if(Dermer_Kachelriess2013_heavy==1)
  {
    cout<<"test Dermer_Kachelriess2013_heavy with energy dispersion"<<endl;

  energy_units="MeV";

  energy_type="kinetic_energy";
  model_name="Dermer_Kachelriess2013_heavy";
  init(energy_type,energy_units,0,model_name);


  E_gamma_band_min_MeV[0]=10;  E_gamma_band_max_MeV[0]=20;
  E_gamma_band_min_MeV[1]=20;  E_gamma_band_max_MeV[1]=30;
  E_gamma_band_min_MeV[2]=30;  E_gamma_band_max_MeV[2]=40;
  E_gamma_band_min_MeV[3]=40;  E_gamma_band_max_MeV[3]=50;
  E_gamma_band_min_MeV[4]=50;  E_gamma_band_max_MeV[4]=60;

  gen_gamma_band_matrices(E_proton_rebin_MeV,
                          E_gamma_band_min_MeV, E_gamma_band_max_MeV,  E_gamma_band_factor);

  name="Shikaze protons kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);
  
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 
  

  name="Shikaze Helium kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);

  helium_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 

  debug=0;// class member

  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum);

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" using precomputed matrices " 
	     <<" model_name="<<model_name
	     <<" energy_dispersion="<<energy_dispersion
	     <<" energy_units="<<energy_units 
	     <<" energy_type ="<<energy_type  
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="      <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
             <<endl;


  valarray<double> emiss_band_nodisp;
  emiss_band_nodisp.resize(emiss_band_.size());
  emiss_band_nodisp=emiss_band_;



  ///

  int energy_dispersion_=1; // above defaults to 1
  init(energy_type,energy_units,0,model_name,energy_dispersion_);



  gen_gamma_band_matrices(E_proton_rebin_MeV,
                          E_gamma_band_min_MeV, E_gamma_band_max_MeV,  E_gamma_band_factor);

  name="Shikaze protons kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);
  
  proton_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 
  

  name="Shikaze Helium kinetic energy MeV";
  cr_spectrum.init(name,parameters, debug);

  helium_spectrum=cr_spectrum.flux(E_proton_rebin_MeV); 

  debug=0;// class member

  emiss_band_=    emissivity_band( proton_spectrum,  helium_spectrum);

  for(int iE_gamma_band=0;  iE_gamma_band<nE_gamma_band; iE_gamma_band++)
  	     cout<<" using precomputed matrices " 
	     <<" model_name="<<model_name
	     <<" energy_dispersion="<<energy_dispersion
	     <<" energy_units="<<energy_units 
	     <<" energy_type ="<<energy_type  
	     <<" Phi="<<Phi<<" MV"
             <<" E_gamma="<< E_gamma_band_min_MeV[iE_gamma_band]<<" -"<< E_gamma_band_max_MeV[iE_gamma_band]
             <<" predicted emiss="      <<emiss_band_[iE_gamma_band] 
             <<" measured = "           <<emiss_measured[iE_gamma_band]
             <<" predicted/measured = " <<emiss_band_[iE_gamma_band]/emiss_measured[iE_gamma_band]
	     <<" edisp/noedisp ="       <<emiss_band_[iE_gamma_band] /emiss_band_nodisp[iE_gamma_band] 
             <<endl;















 }// if  Dermer_Kachelriess2013_heavy==1

  return 0;
};



