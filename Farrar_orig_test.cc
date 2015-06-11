

// development/test version
// A. Strong, MPE
// June 2 2014

using namespace std;
#include<iostream>
#include<cmath>

  ////////////////////////////////////////////////////////

void JanssonFarrarB(double x, double y, double z,
                    double &Breg, double &Bregx, double &Bregy, double &Bregz,
                    double &Bran, double &Branx, double &Brany, double &Branz,
                    int debug )
{
// B-field model of Jansson and  Farrar (2012) ApJ 757, 14 (JF12)
// This implementatin is based on http://www.ccppastroparticle.com/projects/jf12/
// only the rregular field is defined; the random field is included in the parameter list for compatitility with other models, but is set to zero.

  double CGS_U_pi      = acos(-1.);
  double CGS_U_kpc     = 3.08e21;
  double CGS_U_muGauss = 1.e-6;

  double coords_x,coords_y,coords_z;
  // JF12 coordinates have Sun on -x axis. Hence x,y reversed with respect to GALPROP system (= rotation by 180 deg)

  coords_x = -x*CGS_U_kpc ; // uses cm internally 
  coords_y = -y*CGS_U_kpc ; 
  coords_z =  z*CGS_U_kpc ; 


// disk parameters
// values from Table 1 of JF12
double b_arm_1 =  0.1; // micro Gauss field strength in spiral arm region 1 at r=5 kpc
double b_arm_2 =  3.0;
double b_arm_3 = -0.9;
double b_arm_4 = -0.8;
double b_arm_5 = -2.0;
double b_arm_6 = -4.2;
double b_arm_7 =  0.0;
double b_ring  =  0.1; // field strength in molecular ring

double h_disk  = 0.40; // height of transition between disk and toroidal halo
double w_disk  = 0.27; // transition width between disk and toroidal halo

// toroidal halo parameters
double Bn =  1.4 ; // field strength in the north
double Bs = -1.1 ; // ...            in the south
double rn =  9.22; // transition radius in the north
double rs = 16.7 ; // ...               in the south
double wh =  0.20; // transition width
double z0 =  5.3 ; // scale height in z

// X-field parameters
double B0_X   =  4.6;
double Xtheta = 49.0;
double rpc_X  =  4.8; // called r_X^c in paper
double r0_X   =  2.9; // called r_X in paper

// define fixed parameters
double Rmax   = 20*CGS_U_kpc; //    outer boundary of GMF
double rho_GC = 1.*CGS_U_kpc; // interior boundary of GMF

// fixed disk parameters
double inc    = 11.5; // inclination, in degrees
double rmin   = 5.*CGS_U_kpc; // outer boundary of the molecular ring region
double rcent  = 3.*CGS_U_kpc; // inner boundary (disk field is zero within this radius)
double f[8]   = {0.130, 0.165, 0.094, 0.122, 0.13, 0.118, 0.084, 0.156}; // fractions of circumference spanned by each spiral, sums to unity
double rc_B[8] = {5.1, 6.3, 7.1, 8.3, 9.8, 11.4, 12.7, 15.5}; // the radii where the spiral arm boundaries cross the negative x-axis

// x,y,z is a Galactocentric cartesian system, with the Sun on the negative x-axis
double r   = sqrt(coords_x*coords_x + coords_y*coords_y);
double rho = sqrt(coords_x*coords_x + coords_y*coords_y + coords_z*coords_z);
double PHI = atan2(coords_y,coords_x);


  Breg=Bregx=Bregy=Bregz=0.;
  Bran=Branx=Brany=Branz=0.;


// define boundaries outside of which B is zero

if (r > Rmax)     { return ;}
if (rho < rho_GC) { return ;} 


// Disk component: a divergenceless form of Brown et al. (2007)
// 8 spiral regions, 7 free parameters, the 8th set to conserve flux
// set B0 to 1 muG at r=5 kpc 
double B0 = rmin/r*CGS_U_muGauss; 
// the logistic equation, to be multiplied to the toroidal halo field, 
// and (1-zprofile) multiplied to the disk: 
double zprofile      = 1./(1+exp(-2./w_disk*(abs(coords_z)/CGS_U_kpc-h_disk)) ); 
double B_cyl_disk[3] = {0,0,0}; 


// the disk field in cylindrical coordinates 
if ( r > rcent ) // disk field zero elsewhere
{

 if (r < rmin)
 { // circular field in molecular ring
  B_cyl_disk[1] = B0*b_ring*(1-zprofile);
 }
 else
 {
 // use flux conservation to calculate the field strength in the 8th spiral arm
 double bv_B [8] = {b_arm_1, b_arm_2, b_arm_3, b_arm_4, b_arm_5, b_arm_6, b_arm_7, 0.};
 double b7 = 0.;
 for (int i=0; i<7; i++)
 { 
   b7 -= f[i]*bv_B[i]/f[7]; 
 } // last spiral strength is set by the others, to conserve flux 
 bv_B[7] = b7; // NB was on previous comment line in original

 // iteratively figure out which spiral arm the current coordinates (r, phi) corresponds to. 
 double b_disk = 0; 
 double r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI-CGS_U_pi)); 
 if (r_negx > rc_B[7]*CGS_U_kpc) {r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI+  CGS_U_pi)); }
 if (r_negx > rc_B[7]*CGS_U_kpc) {r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI+3*CGS_U_pi)); }

 int ii=-1; // chosen arm

 for (int i=7; i>=0; i--)
 { 
   if (r_negx < rc_B[i]*CGS_U_kpc) { b_disk = bv_B[i]; ii=i;} 
   if(debug==1) cout<<"i="<<i<<" r_negx="<<r_negx<<" rc_B[i]*CGS_U_kpc = "<<  rc_B[i]*CGS_U_kpc <<"  bv_B[i]="<< bv_B[i]  <<" b_disk="<<b_disk<<endl;
 } 
 // "region 8,7,6,..,2" // the disk field in cylindrical coordinates 
 B_cyl_disk[0] = b_disk*B0*sin(CGS_U_pi/180.*inc)*(1-zprofile); 
 B_cyl_disk[1] = b_disk*B0*cos(CGS_U_pi/180.*inc)*(1-zprofile); 

 if(debug==1) if(ii>-1) cout<<" in arm ii="<<ii<<" r_negx="<<r_negx<<" rc_B[ii]*CGS_U_kpc = "<<  rc_B[ii]*CGS_U_kpc <<"  bv_B[ii]="<< bv_B[ii]  <<" b_disk="<<b_disk<<endl;
 
 } 

} 

// --- Toroidal halo component 
double b1, rh; 
double B_h = 0; 
if ( coords_z >= 0) { // NORTH
b1 = Bn*CGS_U_muGauss;
rh = rn; // transition radius between inner-outer region, (units added later)
}
else if ( coords_z < 0 ){ // SOUTH
b1 = Bs*CGS_U_muGauss;
rh = rs;
}
B_h = b1*(1. - 1./(1.+exp(-2./wh*(r/CGS_U_kpc-rh))))*exp(-(abs(coords_z))/(z0*CGS_U_kpc)); // vertical exponentialfall-off
double B_cyl_h[3]={ 0., B_h*zprofile, 0. };

// --- X-field component
// apply units to input parameters
B0_X   *= CGS_U_muGauss;
r0_X   *= CGS_U_kpc;
rpc_X  *= CGS_U_kpc;
Xtheta *= CGS_U_pi/180.;

double rp_X   = 0.; // the mid-plane radius for the field line that pass through r
double B_X    = 0.;
double r_sign = 1.;
if (coords_z<0){ r_sign = -1.;}

// dividing line between region with constant elevation angle, and the interior:

double rc_X = rpc_X + abs(coords_z)/tan(Xtheta);

if (r<rc_X){ // interior region, with varying elevation angle
 rp_X   = r*rpc_X/rc_X;
 B_X    = B0_X * pow(rpc_X/rc_X ,2.) * exp(-rp_X/r0_X);
 Xtheta = atan( abs(coords_z)/ (r-rp_X) ); // modified elevation angle in interior region
 if (coords_z==0.){Xtheta=CGS_U_pi/2.;} // to avoid some NaN
}
else { // exterior region with constant elevation angle
 rp_X = r - abs(coords_z)/tan(Xtheta);
 B_X  = B0_X * rp_X/r * exp(-rp_X/r0_X);
}

// X-field in cylindrical coordinates
double B_cyl_X[3]={ B_X*cos(Xtheta)*r_sign, 0. , B_X*sin(Xtheta) };

// --- add disk + halo components together -------
double B_cyl[3] = {0,0,0};
B_cyl[0] = B_cyl_disk[0] + B_cyl_h[0] + B_cyl_X[0];
B_cyl[1] = B_cyl_disk[1] + B_cyl_h[1] + B_cyl_X[1];
B_cyl[2] = B_cyl_disk[2] + B_cyl_h[2] + B_cyl_X[2];

// compute x,y,z components
// notation here as for other models

  double theta=atan2(y,x); 
  double B_R  =B_cyl[0];
  double B_phi=B_cyl[1];

  Bregx = B_R*cos(theta) - B_phi*sin(theta);
  Bregy = B_R*sin(theta) + B_phi*cos(theta);
  Bregz = B_cyl[2];

  Breg=sqrt(Bregx*Bregx + Bregy*Bregy +  Bregz*Bregz); 

  Bran=sqrt(Branx*Branx + Brany*Brany +  Branz*Branz); // just for consistency

if(debug==1)
{
cout<<" r="<<       r<<" rmin="<<       rmin<<" rcent="<<rcent<<" B0="<<       B0<<" zprofile="<<zprofile<<endl;

cout<<endl;

cout<<"GALPROP coordinates: x = "<<       x <<" y = "<<       y <<" z= "<<       z<<" kpc"<<endl;
cout<<"JF12 coords: x = "<<coords_x<<" y = "<<coords_y<<" z = "<<coords_z<<" cm"<<endl;

cout<<"Bcyl_disk = "<<B_cyl_disk[0]<<"  "<<B_cyl_disk[1]<<"  "<<B_cyl_disk[2]<<endl;
cout<<"Bcyl_h = "<<B_cyl_h[0]<<"  "<<B_cyl_h[1]<<"  "<<B_cyl_h[2]<<endl;
cout<<"Bcyl_X = "<<B_cyl_X[0]<<"  "<<B_cyl_X[1]<<"  "<<B_cyl_X[2]<<endl;
cout<<"Bcyl   = "<<B_cyl[0]<<"          "<<B_cyl[1]<<"         "<<B_cyl[2]<<endl;
cout<<"theta = "<<theta<<" radians"<< endl;
cout  <<"Breg  = "<<Breg   <<" Bregx  = "<<Bregx     <<"  Bregy = "<<Bregy     <<"  Bregz = "<<Bregz     <<endl;
cout <<"Bran  = "<<Bran <<" Branx  = "<<Branx     <<"  Brany = "<<Brany     <<"  Branz = "<<Branz     <<endl;
}

return;
}

int main()
{

  // input  parameters for B_field_3D_model.cc
  double x,y,z; // kpc 
  int debug;
      debug=0;

  // output  parameters for B_field_3D_model.cc
  double Breg, Bregx, Bregy, Bregz;
  double Bran, Branx, Brany, Branz;


 

  for (z=-10; z<=10; z++)
  for (y=-10; y<=10; y++)
  for (x=-10; x<=10; x++)
  {
    JanssonFarrarB(x,  y,  z, 
                   Breg, Bregx, Bregy, Bregz,
                   Bran, Branx, Brany, Branz,
                   debug );

   cout<<"GALPROP coordinates: x = "<<       x <<" y = "<<       y <<" z = "<<       z<<" kpc"
   <<"  Breg   = "<<Breg    <<"  Bregx  = "<<Bregx     <<"  Bregy = "<<Bregy      <<"  Bregz = "<<Bregz     <<endl;
   cout<<"JF12    coordinates: x = "<<      -x <<" y = "<<      -y <<" z = "<<       z<<" kpc"
   <<"  Breg   = "<<Breg   <<"  Bregx  = "<<Bregx     <<"  Bregy = "<<Bregy       <<"  Bregz = "<<Bregz     <<endl;
  }

      return 0;
}
