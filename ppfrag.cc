// g++ frag.f ppfrag.cc -lgfortran -lgfortranbegin

// for normal usage:
void   ppfrag_spec_ini();
double ppfrag_spec_int (double ,double ,int ,int );

using namespace std;
#include<iostream>

///////////////////////////////////////////////////////////////
// interface to ppfrag fortran routines in frag.f
extern "C" void   spec_ini_();
extern "C" double spec_int_ (double* ,double* ,int* ,int* );
///////////////////////////////////////////////////////////////

void   ppfrag_spec_ini()
{

  // reads data from gamfrag.dat and apfrag.dat
  spec_ini_();
  return;
}
///////////////////////////////////////////////////////////////
double ppfrag_spec_int(double ep,double es,int id, int reac)
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


NB energy is TOTAL per nucleon (see frag.f)

  */


  double result=0.;

  // respect limits and avoid stop in fortran routines

  // pp including Kamae
  if(reac==0 && es<ep && id==0)   result=spec_int_(&ep,&es,&id,&reac);

  // other cases including pbar, only >= 10 GeV
  if(reac >0 && es<ep && ep>=10.) result=spec_int_(&ep,&es,&id,&reac);




  return result;
}
////////////////////////////////////////////////////
// test program
////////////////////////////////////////////////////
int test_ppfrag()
{
  cout<<"ppfrag.cc test program"<<endl;

  ppfrag_spec_ini();
  double ep,es;
  int id,reac;

  ep=100.;// proton energy GeV
  es=  1.;// gamma  energy GeV

  id=0;  // 0=photons 1=antiprotons
  reac=0;// reaction

  double result0,result1,result2,result3;
  cout<<"============ ppfrag: gammas"<<endl;

  for (ep= 1.;ep<100.;ep*=1.5)
  for (es= 1.;es< ep ;es*=1.5)
  {
   reac=0;// pp with Kamae
   result0=ppfrag_spec_int(ep,es,id,reac);
   

  
   reac=1; // pp QCSJET only
   result1=ppfrag_spec_int(ep,es,id,reac);

   reac=2;// He-p QCSJET only
   result2=ppfrag_spec_int(ep,es,id,reac);

   reac=3;// p-He QCSJET only
   result3=ppfrag_spec_int(ep,es,id,reac);
   

   cout <<" ep="<<ep<<" gamma es="<<es
        <<"  pp="  <<result0<<"  p-p =" <<result1
        <<"  p-He="<<result2<<" He-p =" <<result3
       <<endl;
  }

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
   

   cout <<" ep="<<ep<<"  pbar es="<<es
        <<"  pp="  <<result0<<"  p-p =" <<result1
        <<"  p-He="<<result2<<" He-p =" <<result3
       <<endl;
  }


  return 0;
}
////////////////////////////////////////////////////
// test driver, enable to use it

int main()
{
  test_ppfrag();
  return 0;
}

