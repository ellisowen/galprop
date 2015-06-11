using namespace std;
#include<vector>
#include <valarray>
#include <string>

class Bremsstrahlung
{
public:
  void   init(string name, valarray<double> parameters, int debug);
  double dSdE(double Ee, double Egamma, int IZ1, int Ne1 );

  void   bremss_emiss(valarray<double> Ee,    valarray<double> Egamma,
                      valarray<double> fluxe,
		      int IZ1, int Ne,
                      valarray<double> &emiss,
                      int method );

  valarray<double>
         bremss_emiss(valarray<double> Ee,    valarray<double> Egamma,
                      valarray<double> fluxe,
		      int IZ1, int Ne,
                      int method );

  double emissivity_band(valarray<double> Ee, valarray<double> fluxe,
                         int IZ1, int Ne,
                         double E_gamma_band_min,
                         double E_gamma_band_max,
                         double E_gamma_factor,
                         int method );

  valarray<double>
         emissivity_band(valarray<double> Ee, valarray<double> fluxe,
                         int IZ1, int Ne,
                         valarray<double> E_gamma_band_min,
                         valarray<double> E_gamma_band_max, 
                         double E_gamma_factor,
                         int method );

  // version using gamma_band matrices precomputed with gen_gamma_band_matrices

  valarray<double>  emissivity_band ( valarray<double> fluxe); 

  int gen_gamma_band_matrices(valarray<double> Ee, 
                              valarray<int> IZ1_list, valarray<int> Ne_list, 
                              valarray<double> target_abundances,
                              valarray<double> E_gamma_band_min,
                              valarray<double> E_gamma_band_max, 
                              double E_gamma_factor,
                              int method );


  void   test(int test_dSdE, int test_emiss);

  string name;
  valarray<double> parameters; 

  double Z,A;
  
private:
  
  int    initialized;
  int    debug;

 // matrices from lepton spectra to gamma-ray bands
  vector<valarray<double> > gamma_band_matrix;
  valarray<double> Ee_grid; // electron energies for matrix
 
};

