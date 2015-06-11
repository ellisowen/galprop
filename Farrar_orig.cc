// copied from http://www.ccppastroparticle.com/projects/jf12/

// disk parameters
double b_arm_1 = ...; // micro Gauss field strength in spiral arm region 1
double b_arm_2 = ...;
double b_arm_3 = ...;
double b_arm_4 = ...;
double b_arm_5 = ...;
double b_arm_6 = ...;
double b_arm_7 = ...;
double b_ring  = ...; // field strength in molecular ring
double h_disk  = ...; // height of transition between disk and toroidal halo
double w_disk  = ...; // transition width between disk and toroidal halo

// toroidal halo parameters
double Bn = ...; // field strength in the north
double Bs = ...; // ... in the south
double rn = ...; // transition radius in the north
double rs = ...; // ... in the south
double wh = ...; // transition width
double z0 = ...; // scale height in z

// X-field parameters
double B0_X   = ...;
double Xtheta = ...;
double rpc_X  = ...; // called r_X^c in paper
double r0_X   = ...; // called r_X in paper

// define fixed parameters
double Rmax   = 20*CGS_U_kpc; // outer boundary of GMF
double rho_GC = 1.*CGS_U_kpc; // interior boundary of GMF

// fixed disk parameters
double inc    = 11.5; // inclination, in degrees
double rmin   = 5.*CGS_U_kpc; // outer boundary of the molecular ring region
double rcent  = 3.*CGS_U_kpc; // inner boundary (disk field is zero within this radius)
double f[8]   = {0.130, 0.165, 0.094, 0.122, 0.13, 0.118, 0.084, 0.156}; // fractions of circumference spanned by each spiral, sums to unity
double rc_B[8] = {5.1, 6.3, 7.1, 8.3, 9.8, 11.4, 12.7, 15.5}; // the radii where the spiral arm boundaries cross the negative x-axis

// x,y,z is a Galactocentric cartesian system, with the Sun on the negative x-axis
double r   = sqrt(coords.x*coords.x + coords.y*coords.y);
double rho = sqrt(coords.x*coords.x + coords.y*coords.y + coords.z*coords.z);
double PHI = atan2(coords.y,coords.x);
double z   = coords.z;

// define boundaries outside of which B is zero
if (r > Rmax) { return vec3(0,0,0);}
if (rho < rho_GC) { return vec3(0,0,0);} 

// Disk component: a divergenceless form of Brown et al. (2007)
// 8 spiral regions, 7 free parameters, the 8th set to conserve flux
// set B0 to 1 muG at r=5 kpc 
double B0 = rmin/r*CGS_U_muGauss; 
// the logistic equation, to be multiplied to the toroidal halo field, 
// and (1-zprofile) multiplied to the disk: 
double zprofile      = 1./(1+exp(-2./w_disk*(abs(z)/CGS_U_kpc-h_disk)) ); 
double B_cyl_disk[3] = {0,0,0}; 


// the disk field in cylindrical coordinates 
if ( r > rcent ) // disk field zero elsewhere
{
 if (r < rmin) { // circular field in molecular ring
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
 } // last spiral strength is set by the others, to conserve flux bv_B[7] = b7; 

 // iteratively figure out which spiral arm the current coordinates (r, phi) corresponds to. 
 double b_disk = 0; 
 double r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI-CGS_U_pi)); 
 if (r_negx > rc_B[7]*CGS_U_kpc) {r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI+CGS_U_pi)); }
 if (r_negx > rc_B[7]*CGS_U_kpc) {r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI+3*CGS_U_pi)); }
 for (int i=7; i>=0; i--){ 
    if (r_negx < rc_B[i]*CGS_U_kpc) { b_disk = bv_B[i];} } 
 // "region 8,7,6,..,2" // the disk field in cylindrical coordinates 
 B_cyl_disk[0] = b_disk*B0*sin(CGS_U_pi/180.*inc)*(1-zprofile); 
 B_cyl_disk[1] = b_disk*B0*cos(CGS_U_pi/180.*inc)*(1-zprofile); 
 } 
} 

// --- Toroidal halo component 
double b1, rh; 
double B_h = 0; 
if ( z >= 0) { // NORTH
b1 = Bn*CGS_U_muGauss;
rh = rn; // transition radius between inner-outer region, (units added later)
}
else if ( z < 0 ){ // SOUTH
b1 = Bs*CGS_U_muGauss;
rh = rs;
}
B_h = b1*(1. - 1./(1.+exp(-2./wh*(r/CGS_U_kpc-rh))))*exp(-(abs(z))/(z0*CGS_U_kpc)); // vertical exponentialfall-off
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
if (z<0){ r_sign = -1.;}

/ dividing line between region with constant elevation angle, and the interior:
double rc_X = rpc_X + abs(z)/tan(Xtheta);

if (r<rc_X){ // interior region, with varying elevation angle
 rp_X   = r*rpc_X/rc_X;
 B_X    = B0_X * pow(rpc_X/rc_X ,2.) * exp(-rp_X/r0_X);
 Xtheta = atan( abs(z)/ (r-rp_X) ); // modified elevation angle in interior region
 if (z==0.){Xtheta=CGS_U_pi/2.;} // to avoid some NaN
}
else { // exterior region with constant elevation angle
 rp_X = r - abs(z)/tan(Xtheta);
 B_X  = B0_X * rp_X/r * exp(-rp_X/r0_X);
}

// X-field in cylindrical coordinates
double B_cyl_X[3]={ B_X*cos(Xtheta)*r_sign, 0. , B_X*sin(Xtheta) };

// --- add disk + halo components together -------
double B_cyl[3] = {0,0,0};
B_cyl[0] = B_cyl_disk[0] + B_cyl_h[0] + B_cyl_X[0];
B_cyl[1] = B_cyl_disk[1] + B_cyl_h[1] + B_cyl_X[1];
B_cyl[2] = B_cyl_disk[2] + B_cyl_h[2] + B_cyl_X[2];

