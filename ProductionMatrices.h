
//////////////////////////////////////////////////////
class ProductionMatrices
{
public:
  int nE_proton;
  int nE_gamma;

  valarray<double> E_proton, E_gamma;


  // original Huang matrices

  valarray<double>          prodxsection_p_matrix; // total inelastic cross-section, mb
  valarray<double>          prodxsection_he_matrix;// total inelastic cross-section, mb

  vector<valarray<double> > gamma_decay_p_matrix;  // differential gamma-ray spectrum, GeV-1
  vector<valarray<double> > gamma_decay_he_matrix; // differential gamma-ray spectrum, GeV-1


  // original FLUKA matrices
  int nE_proton_FLUKA;
  int nE_gamma_FLUKA;

  valarray<double> E_proton_FLUKA, E_gamma_FLUKA;
  vector<valarray<double> > gamma_p_FLUKA_matrix;   // differential gamma-ray cross-section, mb GeV-1
  vector<valarray<double> > gamma_he_FLUKA_matrix;  // differential gamma-ray cross-section, mb GeV-1 NOT YET USED


  // rebinned matrices for differential spectra

  int nE_proton_rebin;
  int nE_gamma_rebin;
  valarray<double> E_proton_rebin, E_gamma_rebin;
  vector<valarray<double> > gamma_decay_p_matrix_rebin;
  vector<valarray<double> > gamma_decay_he_matrix_rebin;

  valarray<double>          gamma_decay_p_emissivity_rebin;
  valarray<double>          gamma_decay_he_emissivity_rebin;

  // matrices from p, He spectra to gamma-ray bands
  vector<valarray<double> > gamma_band_p_matrix;
  vector<valarray<double> > gamma_band_he_matrix;

  string energy_units; // "GeV" or "MeV"  (per nucleon)
  string energy_type;  // "total_energy" or "kinetic_energy" (per nucleon)
  string model_name;   // "Huang" or "Kachelriess"
  int    energy_dispersion; // 1 = use energy dispersion
  
  int init(string energy_type, string energy_units, int debug=0, string model_name="Huang", int energy_dispersion=0);

  int read(string input_directory=".");
  int print();
  int rebin(valarray<double> E_proton_rebin, valarray<double> E_gamma_rebin);
  int rebin_with_energy_dispersion(valarray<double> E_proton_rebin, valarray<double> E_gamma_rebin, valarray<double> E_gamma_true  );

  int emissivity(valarray<double> proton_spectrum_rebin, valarray<double> helium_spectrum_rebin );

  valarray<double>  emissivity_total(valarray<double> proton_spectrum_rebin, valarray<double> helium_spectrum_rebin );


  double emissivity_band(valarray<double> proton_spectrum_rebin, valarray<double> helium_spectrum_rebin,
                         valarray<double> E_proton_rebin,
                         double E_gamma_band_min, double E_gamma_band_max, double E_gamma_factor);


  valarray<double>  emissivity_band
          (valarray<double> proton_spectrum_rebin, valarray<double> helium_spectrum_rebin,
           valarray<double> E_proton_rebin_,
           valarray<double>E_gamma_band_min, valarray<double> E_gamma_band_max, double E_gamma_factor);

  int gen_gamma_band_matrices(valarray<double> E_proton_rebin,
                              valarray<double>E_gamma_band_min, valarray<double> E_gamma_band_max, double E_gamma_factor);

  // version using gamma_band matrices precomputed with gen_gamma_band_matrices
  valarray<double>  emissivity_band
                   (valarray<double> proton_spectrum_rebin, valarray<double> helium_spectrum_rebin);                         

  // Kachelriess & Ostapchenko 2012 functions
  void   ppfrag_spec_ini();
  double ppfrag_spec_int(double ep,double es,int id, int reac);

  // Kachelriess & Ostapchenko 2013 functions
  double ppfrag_v02     (double ep,double es,int id, int reac, int kamae_QGS); //AWS20130105

  int test_ppfrag();




  int test();

 private:
  int debug;
  int matrices_read;
  int matrices_rebinned;

  int ppfrag_read;

  EnergyDispersion energyDispersion;

};
