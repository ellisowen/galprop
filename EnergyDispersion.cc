using namespace std;
#include<iostream>
#include<fstream>
#include<cstring>
#include<cstdlib>
#include<vector>
#include<valarray>

#include"fitsio.h"

#include"EnergyDispersion.h"

//////////////////////////////////////////////////////////////////

void EnergyDispersion::read(int debug)
{
 

  cout<<">> EnergyDispersion::read"<<endl;

  fitsfile *fptr;
  int status=-1;
  char *Fermi_rsp_file=new char[100];
 
  strcpy (Fermi_rsp_file,"rsp_big.fits");
  cout<<"EmissivityData:: ReadEnergyDispersion"<<endl;
  cout<<"Fermi rsp file: "<<Fermi_rsp_file<<endl;

  fits_open_file(&fptr,Fermi_rsp_file,READONLY,&status);
  if(fptr==NULL||status!=0)
   {
      cout<<"no Fermi rsp file called "<<Fermi_rsp_file; return ;
   }

  fits_movnam_hdu(fptr,ANY_HDU,"MATRIX",0, &status);
  cout<<"Fermi rsp file: HDU MATRIX status="<<status<<endl;


  double *fitsdata; 
  long nrows=1;
  int ncolumns=1;
  fits_get_num_rows(fptr,&nrows,&status);
  cout<<"Fermi rsp file HDU MATRIX nrows="<<nrows<<endl;

  fits_get_num_cols(fptr,&ncolumns,&status);
  cout<<"Fermi rsp file HDU MATRIX top level ncolumns="<<ncolumns<<endl;

  int colnum=0;
  fits_get_colnum(fptr,CASESEN,"MATRIX",&colnum,&status);
  cout<<"Fermi rsp file MATRIX colnum="<<colnum<<endl;

  int  typecode=0;
  long repeat=0;
  long width =0;
  fits_get_coltype(fptr,colnum, &typecode,&repeat,&width,&status);
  cout<<" Fermi rsp file MATRIX fits_get_coltype     typecode="<<typecode<<endl;
  cout<<" Fermi rsp file MATRIX fits_get_coltype       repeat="<<repeat  <<endl;
  cout<<" Fermi rsp file MATRIX fits_get_coltype       width ="<<width   <<endl;
  cout<<" Fermi rsp file MATRIX fits_get_coltype       status="<<status  <<endl;

  // variable length rows require this routine
  long offset=0;
  ncolumns=0;
  for (LONGLONG irow=1;irow<=nrows;irow++)
  {
   fits_read_descript(fptr,colnum,irow, &repeat,&offset,&status);

   if(debug==1)
   {
    cout<< "irow="<<irow;
    cout<<" Fermi rsp file MATRIX  fits_read_descript       offset="<<offset  <<endl;
    cout<<" Fermi rsp file MATRIX  fits_read_descript       repeat="<<repeat  <<endl;
    cout<<" Fermi rsp file MATRIX  fits_read_descript       status="<<status  <<endl;
   }

   if(repeat>ncolumns)ncolumns=repeat; // get max length this way
  }

  cout<<"Fermi rsp file MATRIX ncolumns from max repeat="<<ncolumns<<endl;
  

  fitsdata=new double[nrows];
  LONGLONG nelements=nrows;
  LONGLONG first_row    =1;
  LONGLONG first_element=1;

  EnergyDispersionMatrix.resize(nrows);
  for(int i=0;i<nrows;i++)
  {
   EnergyDispersionMatrix[i].resize(nrows);
   for(int j=0;j<ncolumns;j++)
   {   
              EnergyDispersionMatrix[i][j]=0;
   }
  }



  int i,j;



  for (int irow=1;irow<=nrows;irow++)
  {

  double nulval=0;
  int    anynul=0;
  for(int k=0;k<nrows;k++) fitsdata[k]=0; // since may be undefined values

  fits_read_descript(fptr,colnum,irow, &repeat,&offset,&status);
  nelements=repeat;
  first_row=irow;

  fits_read_col(fptr,TDOUBLE,colnum,first_row,first_element,nelements,&nulval,fitsdata,&anynul,&status);


  if(debug==1)
  {
  cout<<"Fermi rsp file: MATRIX read status="<<status<<endl;
  cout<<"Fermi rsp file: MATRIX read anynul="<<anynul<<endl;


   cout<<"fitsdata irow= "<<irow<<" nelements= "<<nelements <<":";
   for(int k=0;k<nelements;k++)                                              cout<<fitsdata[k]<<" ";
   cout<<endl;
  }

  for(int k=0;k<nelements;k++){i=irow-1;j=k;         EnergyDispersionMatrix[i][j]=fitsdata[k];}


  } // irow


 if(debug==1)
 {
  cout<<"EnergyDispersionMatrix: filled to i="<<i<<",  j="<<j<<endl;


  for (int i=0;i<nrows;i++)
  {
    cout<<"EnergyDispersionMatrix as read: i="<<i<<": ";
   for(int j=0;j<ncolumns;j++)
   {
     cout<<  EnergyDispersionMatrix[i][j]  <<" ";
 
   }
   cout<<endl;
  }
}// debug
  

 

 // normalize each row. NB assumes the normalization is 1, i.e. the dispersion matrix is complete
  for (int i=0;i<nrows;i++) EnergyDispersionMatrix[i] /= EnergyDispersionMatrix[i].sum();


  // notation used in class
  nE_true=nrows;
  nE_meas=ncolumns;
 
 if(debug==1)
 {
  for (int i=0;i<nE_true;i++)
  {
    cout<<"normalized  EnergyDispersionMatrix: i="<<i<<": ";
   for(int j=0;j<nE_meas;j++)
   {
     cout<<  EnergyDispersionMatrix[i][j]  <<" ";

 
   }
    cout<< " sum="<< EnergyDispersionMatrix[i].sum();
    cout<<endl;
  }
 }// debug


 // following must be done using EBOUNDS and MATRIX
 E_true.resize(nE_true);
 E_meas.resize(nE_meas);

 E_true[0]        =20.;       // MeV until read from dataset
 E_true[nE_true-1]=2.029176e5;// MeV until read from dataset
 E_meas[0]        =20.;       // MeV until read from dataset
 E_meas[nE_meas-1]=2.029176e5;// MeV until read from dataset




 log_E_true_min=log(E_true[0]);
 log_E_meas_min=log(E_meas[0]);
 dlog_E_true=log( E_true[nE_true-1]/ E_true[0])/(nE_true-1);
 dlog_E_meas=log( E_meas[nE_meas-1]/ E_meas[0])/(nE_meas-1);

 for (int i=0;i<nE_true;i++) E_true[i]=exp(log_E_true_min+i* dlog_E_true);
 for (int j=0;j<nE_meas;j++) E_meas[j]=exp(log_E_meas_min+j* dlog_E_meas);

 dE_true.resize(nE_true);
 dE_meas.resize(nE_meas);
                            dE_true[0]=E_true[1]  -E_true[0];
 for (int i=1;i<nE_true;i++)dE_true[i]=E_true[i]  -E_true[i-1];

                            dE_meas[0]=E_meas[1]  -E_meas[0];
 for (int j=1;j<nE_meas;j++)dE_meas[j]=E_meas[j]  -E_meas[j-1];



 cout<<   "EnergyDispersion read complete, nE_true="<<nE_true<<" nE_meas="<<nE_meas<<endl;
 cout<<   "EnergyDispersion read complete,  log_E_true_min="<<log_E_true_min<<" log_E_meas_min="<<log_E_meas_min<<endl;
 cout<<   "EnergyDispersion read complete, dlog_E_true="<<dlog_E_true<<" dlog_E_meas="<<dlog_E_meas<<endl;
 
if(debug==1)
 {
  for (int i=1;i<nE_true;i++)cout<<"i="<<i<<" E_true[i]="<< E_true[i]<<" dE_true[i]="<< dE_true[i]<<endl;
  for (int j=1;j<nE_meas;j++)cout<<"j="<<j<<" E_meas[j]="<< E_meas[j]<<" dE_meas[j]="<< dE_meas[j]<<endl;
 }



 for (int i=0;i<nrows;i++) EnergyDispersionMatrix[i] /= dE_meas; // units MeV^-1




 if(debug==1)
 {
  for (int i=0;i<nE_true;i++)
  {
    cout<<"  EnergyDispersionMatrix MeV^-1: i="<<i<<": ";
    double sum=0;
   for(int j=0;j<nE_meas;j++)
   {
     cout<<  EnergyDispersionMatrix[i][j]  <<" ";
     sum +=  EnergyDispersionMatrix[i][j] * dE_meas[j];
 
   }
   cout<< " sum with MeV ="<< sum;
    cout<<endl;
  }
 }// debug


 initialized=1;

 cout<<"<< EnergyDispersion::read"<<endl;

  return;
};

////////////////////////////////////////////////////////////////////////////

double EnergyDispersion::value(double E_true,double E_meas, int debug)
{ 
  if(initialized!=1){cout<<" not initialized: initializing"<<endl; read(debug);}

  int iE_true=(log(E_true)-log_E_true_min)/dlog_E_true;
  int iE_meas=(log(E_meas)-log_E_meas_min)/dlog_E_meas;

  if(debug==1)cout<<"  E_true="<< E_true<<"  E_meas="<< E_meas;
  if(debug==1)cout<<" iE_true="<<iE_true<<" iE_meas="<<iE_meas;

  double value_=0.;

  if(iE_true<0||iE_true>nE_true-1  
   ||iE_meas<0||iE_meas>nE_meas-1) 
  {
   if(debug==1)cout<<"  EnergyDispersionMatrix outside range,  value=" << value_ <<endl;
   return value_;
  }

  value_ = EnergyDispersionMatrix[iE_true][iE_meas];

  if(debug==1)cout<<"  EnergyDispersionMatrix value=" << value_ <<endl;

  return value_;
}
////////////////////////////////////////////////////////////////////////////
  void  EnergyDispersion::ApplyEnergyDispersion(valarray<double> E,  valarray<double> &spectrum, int debug)
{
  // in-place application of energy dispersion to spectrum
  int nE=spectrum.size();
  valarray<double> work;
  work.resize(nE);  
  work = 0.;


  for (int j=0;j<nE;j++)
  for (int i=0;i<nE;i++)
  {
     double E_true_=E[i];
     double E_meas_=E[j];
     double value_=value(E_true_,E_meas_,debug); // matrix element per MeV of E_meas
     work[j] += spectrum[i] * value_ * E_true_;  // E_true_: since log integration

     //    cout<<"ApplyEnergyDispersion : i="<<i<<" E_true_="<<E_true_<<" j="<<j<<"  E_meas_="<<E_meas_<< " spectrum input = "<<spectrum[i]<<" value_="<<value_<<" output =  "<<work[j]<<endl;
  }

  work *= log(E[1]/E[0]); // log integration

  double  input_sum=0;
  double output_sum=0;
  for (int i=0;i<nE;i++){input_sum+=spectrum[i]*E[i];  output_sum+=work[i]*E[i];}
   input_sum *= log(E[1]/E[0]); // log integration
  output_sum *= log(E[1]/E[0]); // log integration

  for (int i=0;i<nE;i++)cout<<"ApplyEnergyDispersion : i="<<i<<"  E[i]="<<E[i]<<" spectrum input = "<<spectrum[i] <<" output = "<<work[i]
                            <<" integration input="<<input_sum<<" output="<<output_sum<<" output/input integrals="<<output_sum/input_sum<<endl;
  
  spectrum = work; // in-place result

  return;
}
////////////////////////////////////////////////////////////////////////////
  void  EnergyDispersion::ApplyEnergyDispersion(valarray<double> E_true_,  valarray<double>  spectrum_true,valarray<double> E_meas_,  valarray<double> &spectrum_meas, int debug)
{
  //  application of energy dispersion to spectrum with independent E_true_, E_meas_
  // NB watch for confusion with E_true,E_meas used in the matrix !

  int nE_true=spectrum_true.size();
  int nE_meas=spectrum_meas.size();

    
  spectrum_meas = 0.;


  for (int j=0;j<nE_meas;j++)
  for (int i=0;i<nE_true;i++)
  {
     
     double value_     = value(E_true_[i], E_meas_[j], debug);    // matrix element per MeV of E_meas_
     spectrum_meas[j] += spectrum_true[i] * value_ * E_true_[i];  // E_true_: since log integration

     if(debug==1)      cout<<"ApplyEnergyDispersion : i="<<i<<" E_true_="<<E_true_[i]<<" j="<<j<<"  E_meas_="<<E_meas_[j]<< " spectrum input = "<<spectrum_true[i]<<" value_="<<value_<<" output =  "<<spectrum_meas[j]<<endl;
  }

    spectrum_meas *= log(E_true_[1]/E_true_[0]); // log integration

 

  double  input_sum=0;
  double output_sum=0;

  for (int i=0;i<nE_true;i++) input_sum+=spectrum_true[i]*E_true_[i];
  for (int j=0;j<nE_meas;j++)output_sum+=spectrum_meas[j]*E_meas_[j];  

   input_sum *= log(E_true_[1]/E_true_[0]); // log integration
  output_sum *= log(E_meas_[1]/E_meas_[0]); // log integration


 if(debug==1)
 {
  for (int i=0;i<nE_true;i++)cout<<"ApplyEnergyDispersion : i="<<i<<"  E_true_[i]="<<E_true_[i]<<" spectrum input  = "<<spectrum_true[i] <<endl;
  for (int j=0;j<nE_meas;j++)cout<<"ApplyEnergyDispersion : j="<<j<<"  E_meas_[j]="<<E_meas_[j]<<" spectrum output = "<<spectrum_meas[j] <<endl;
    
  cout<<"ApplyEnergyDispersion :  integration input="<<input_sum<<" output="<<output_sum<<" output/input integrals="<<output_sum/input_sum<<endl;
  cout<<" true sum="<<spectrum_true.sum()<<" meas sum="<<spectrum_meas.sum()<<" dlogtrue="<<  log(E_true_[1]/E_true_[0])<<" dlogmeas="<<  log(E_meas_[1]/E_meas_[0])  <<endl;
 }

  return;
}
////////////////////////////////////////////////////////////////////////////
void EnergyDispersion::test()
{
  cout<<">>  EnergyDispersion::test"<<endl;

  EnergyDispersion energyDispersion;
  int debug=0;
  energyDispersion.read(debug); // remove to test initialization check

  double E_true=100;
  double E_meas=200;
  double E_meas_factor=1.050;
  debug=0;
  double value;
  for(E_true=10;E_true<1e6;E_true*=2)
  {
   double sum=0;
   double sum1=0;
  for(E_meas=10;E_meas<1e6;E_meas*=E_meas_factor)
  {
   value=energyDispersion.value(E_true,E_meas,debug);
   double dE_meas=E_meas*(sqrt(E_meas_factor) - 1./sqrt(E_meas_factor) );
   sum +=value * dE_meas;
   sum1+=value *  E_meas;
   cout<<"  E_true="<< E_true<<"  E_meas="<< E_meas<<" value="<<value<<endl;
  }
  sum1*=log(E_meas_factor); // log integration more accurate
  cout<<"E_true="<<E_true<<"  sum (value*dE_meas) ="<<sum<<"    log (E factor) * sum (value*E_meas) = "<<sum1<<endl;

 }

  int test_in_place=0;
  if(test_in_place==1)
  {
  // test in-place dispersion
  int nE=300;
  valarray<double> E,spectrum;
  E.resize(nE);
  spectrum.resize(nE);
  for(int i=0;i<nE;i++) E[i]=10*pow(10,i*.01);
 
  debug=1;


  spectrum=0;
  spectrum[ 10]=1;
  spectrum[ 20]=1;
  spectrum[ 30]=1;
  spectrum[ 40]=1;
  spectrum[ 50]=1;
  spectrum[ 60]=1;
  spectrum[ 70]=1;
  spectrum[100]=1;
  spectrum[150]=1;
  spectrum[200]=1;

  ApplyEnergyDispersion(E,spectrum,debug);

  spectrum=1;
  ApplyEnergyDispersion(E,spectrum,debug);

  // shows loss at low energies for steep power law
 for(int i=0;i<nE;i++) spectrum[i]=pow(E[i],-3);
 ApplyEnergyDispersion(E,spectrum,debug);

 // start at higher energies, less loss
 for(int i=0;i<nE;i++) E[i]=200*pow(10,i*.01);
 for(int i=0;i<nE;i++) spectrum[i]=pow(E[i],-3);
 ApplyEnergyDispersion(E,spectrum,debug);


// wide energies but spectrum starts higher so covered. integrals agree to 1%
 for(int i=0;i<nE;i++) E[i]=10*pow(10,i*.01);
 spectrum=0;
 for(int i=0;i<nE;i++)if(E[i]>50.) spectrum[i]=pow(E[i],-3);
 ApplyEnergyDispersion(E,spectrum,debug);

  }// test in-place


 // test general dispersion
 // case where energies equal
  debug=1;
  int nE_true=200;
  int nE_meas=200;
  valarray<double> E_true_,E_meas_,spectrum_true,spectrum_meas;
  E_true_.resize(nE_true);
  E_meas_.resize(nE_meas);
  spectrum_true.resize(nE_true);
  spectrum_meas.resize(nE_meas);
  for(int i=0;i<nE_true;i++) E_true_[i]=10*pow(10,i*.01);
  for(int j=0;j<nE_meas;j++) E_meas_[j]=10*pow(10,j*.01);
  spectrum_true=0;
  spectrum_true[100]=1;
  ApplyEnergyDispersion(E_true_,spectrum_true,E_meas_,spectrum_meas,debug);

  
  // case where Emeas same start but different bins
  for(int i=0;i<nE_true;i++) E_true_[i]= 10*pow(10,i*.01 );
  for(int j=0;j<nE_meas;j++) E_meas_[j]= 10*pow(10,j*.007);
  spectrum_true=0;
  spectrum_true[100]=1;
  ApplyEnergyDispersion(E_true_,spectrum_true,E_meas_,spectrum_meas,debug);

  
  // case where Emeas same start but different bins
  for(int i=0;i<nE_true;i++) E_true_[i]= 10*pow(10,i*.01 );
  for(int j=0;j<nE_meas;j++) E_meas_[j]= 10*pow(10,j*.020);
  spectrum_true=0;
  spectrum_true[100]=1;
  ApplyEnergyDispersion(E_true_,spectrum_true,E_meas_,spectrum_meas,debug);

    // case constant spectrum, measured covers more
  for(int i=0;i<nE_true;i++) E_true_[i]= 50*pow(10,i*.01 );
  for(int j=0;j<nE_meas;j++) E_meas_[j]= 50*pow(10,j*.020);
  spectrum_true=1;
  
  ApplyEnergyDispersion(E_true_,spectrum_true,E_meas_,spectrum_meas,debug);

  // steep power law covered by measured, might show loss due to response since starts at 20 MeV
  for(int i=0;i<nE_true;i++) E_true_[i]= 10*pow(10,i*.01 );
  for(int j=0;j<nE_meas;j++) E_meas_[j]= 10*pow(10,j*.020);
  spectrum_true=0;
  for(int i=0;i<nE_true;i++)if(E_true_[i]>20.) spectrum_true[i]=pow(E_true_[i],-3);
  ApplyEnergyDispersion(E_true_,spectrum_true,E_meas_,spectrum_meas,debug);

  // steep power law covered by measured, shows no  loss due to response since starts at 300 MeV
  for(int i=0;i<nE_true;i++) E_true_[i]= 10*pow(10,i*.01 );
  for(int j=0;j<nE_meas;j++) E_meas_[j]= 10*pow(10,j*.020);
  spectrum_true=0;
  for(int i=0;i<nE_true;i++)if(E_true_[i]>300.) spectrum_true[i]=pow(E_true_[i],-3);
  ApplyEnergyDispersion(E_true_,spectrum_true,E_meas_,spectrum_meas,debug);



 cout<<"<<  EnergyDispersion::test"<<endl;
 return;
}

///////////////////////////////////main test program, normally commented out ////////////////////////////////////////
// g++ EnergyDispersion.cc -I/afs/ipp-garching.mpg.de/home/a/aws/propagate/c/cfitsio/3.26/gcc_sles11_olga2/cfitsio/include -L/afs/ipp-garching.mpg.de/home/a/aws/propagate/c/cfitsio/3.26/gcc_sles11_olga2/cfitsio/lib -lcfitsio -o energydispersion

/*
int main()
{
  EnergyDispersion energyDispersion;
  energyDispersion.test();
  return 0;
}
*/
