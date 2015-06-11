using namespace std;
#include <cmath>
#include <iostream>
#include "constants.h"

// hadronic energy losses for nuclei
// AWS20140205
// first draft as placeholder, using simple formulation, while more accurate formulation is under study
// only valid for p and He so far, and only a rough approximation
// only p target so far


///////////////////////////////////////////////

double hadronic_energy_loss(int Z, int A, double Ekin, double nhcm3, double he_to_h)
  
{
  double result = 0.;

  double Ekin_threshold= 400.; //threshold for pion production, MeV/nucleon, roughly

  if(Ekin<Ekin_threshold) return result;

  // 44 mb A^0.7
  double total_inelastic_cross_section = 44.e-27 * pow(A,0.7); // cm^2

  double inelasticity = 0.6;

  // dE/dt in MeV s^-1  units :  E * cm^2 * cm^-3 * cm s^-1

  result = Ekin * total_inelastic_cross_section * inelasticity * nhcm3 * C;

  return result;
}

/////////////////////////////////////////////////////// test program

int main()
{ 
  int Z=1;
  int A=1;
  double Ekin=1e3;
  double nhcm3=1;
  double he_to_h=0.1;

  for (A=1;A<=26;A++)
  for (Ekin=1;Ekin<1e6;Ekin*=10)
  {
   double dEdt= hadronic_energy_loss(Z, A, Ekin, nhcm3, he_to_h);

   cout<<"A="<<A<<" Ekin="<<Ekin<<" dEdt="<<dEdt<<" MeV s^-1"<<endl;
  }


}
