using namespace std;
#include<vector>
#include <valarray>

class InverseCompton
{
public:
  void   init(int debug);
  double IC_sigma_KN          (double Ee, double Eph, double eta, double Egamma);
  double IC_sigma_KN_isotropic(double Ee, double Eph,             double Egamma);

  void   IC_emiss(valarray<double> Ee,   valarray<double> Eph, double eta,valarray<double> Egamma,
                  valarray<double> fluxe,valarray<double> nph,            valarray<double> &emiss,
                  int method );

  void   IC_loss (         double Ee,  
                  valarray<double> Eph,  valarray<double> nph,    
		           double &dEdt, double &mean_energy_loss,
                  int method);
                  

  void   test(int options);

private:
 
  double constant;
  int    initialized;
  int    debug;
};
