class EnergyDispersion
{
 public:
  int nE_true,nE_meas;
  vector<valarray<double> > EnergyDispersionMatrix;
         valarray<double>     E_true,         E_meas;
	          double  log_E_true_min, log_E_meas_min;
	 valarray<double>    dE_true    ,    dE_meas    ;
		  double dlog_E_true,    dlog_E_meas;
  void read(int debug);
  double value(double E_true,double E_meas, int debug);

  void ApplyEnergyDispersion(valarray<double> E,  valarray<double> &spectrum, int debug);
  void ApplyEnergyDispersion(valarray<double> E_true_,  valarray<double>  spectrum_true,valarray<double> E_meas_,  valarray<double> &spectrum_meas, int debug);

  void test();
  int initialized;
};
