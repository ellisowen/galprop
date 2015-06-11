using namespace std;
#include<vector>
#include <valarray>
#include <string>

class CR_Spectrum
{
public:
  void             init(string name, valarray<double> parameters, int debug);
  valarray<double>  flux(valarray<double> E);

  valarray<double> fluxe(valarray<double> Ee, double r);// for stellar case
  
  void   test();

  string name;
  valarray<double> parameters; 

  double Z,A;
  
private:
  
  int    initialized;
  int    debug;
};
