//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * propel.cc *                                  galprop package * 10/18/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// for method see notebook of 970300

#include <iostream>
#include <sstream>

#include "galprop_classes.h"
#include "galprop_internal.h"

using namespace std;

#include <ErrorLogger.h>

int propel_arrays_initialize=1;

Distribution alpha1_r, alpha1_x, alpha1_y, alpha1_z, alpha1_p;
Distribution alpha2_r, alpha2_x, alpha2_y, alpha2_z, alpha2_p;
Distribution alpha3_r, alpha3_x, alpha3_y, alpha3_z, alpha3_p;
Distribution Nr1_, Nx1_, Ny1_, Nz1_, Np1_;
Distribution Nr2_, Nx2_, Ny2_, Nz2_, Np2_;
Distribution Nr3_, Nx3_, Ny3_, Nz3_, Np3_;
Distribution total_source_function;
Distribution Work;

/* may replace with double*
valarray<double> particle_cr_density_arr;                                               //AWS20110126
valarray<double> total_source_function_arr;                                             //AWS20110126
valarray<double> alpha1_r_arr, alpha1_x_arr, alpha1_y_arr, alpha1_z_arr, alpha1_p_arr;  //AWS20110126
valarray<double> alpha2_r_arr, alpha2_x_arr, alpha2_y_arr, alpha2_z_arr, alpha2_p_arr;  //AWS20110126
valarray<double> alpha3_r_arr, alpha3_x_arr, alpha3_y_arr, alpha3_z_arr, alpha3_p_arr;  //AWS20110126
valarray<double> Work_arr;                                                              //AWS20110126
*/

double *particle_cr_density_arr;                                                        //AWS20110128
double *total_source_function_arr;                                                      //AWS20110128
double *alpha1_r_arr,*alpha1_x_arr,*alpha1_y_arr,*alpha1_z_arr,*alpha1_p_arr;           //AWS20110128
double *alpha2_r_arr,*alpha2_x_arr,*alpha2_y_arr,*alpha2_z_arr,*alpha2_p_arr;           //AWS20110128
double *alpha3_r_arr,*alpha3_x_arr,*alpha3_y_arr,*alpha3_z_arr,*alpha3_p_arr;           //AWS20110128
double *Work_arr;                                                                       //AWS20110128

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galprop::propel(Particle& particle) {

  INFO("Entry");

  int start=1;
  double end_timestep_sec = galdef.end_timestep*year2sec;
  double dt = galdef.start_timestep*year2sec;
  double factor = dt*pow(kpc2cm, -2.);  
  double t = 0.;
  int converged = 0;   //AWS20110114
  int method;          //AWS20111018  1=Crank-Nicolson, 2=fully time explicit 
  int turbo;           //AWS20110127  0= use Distributions, 1=use linear arrays (method 2 only)
       

  if(galdef.solution_method!=1 && galdef.solution_method!=2 && galdef.solution_method!=21   ) return 0; //AWS20110127
  
  
  propel_diagnostics(); // AWS20110113 initialization for new species

  ostringstream buf;

  if (0 == particle.primary_source_function.max() && 
      0 == particle.secondary_source_function.max()) {

    buf.str("");
    buf << "Zero primary and secondary source functions. No propagation will be done";
    WARNING(buf.str());
    INFO("Exit");
    return 0;
 
  }

  if (2 == particle.n_spatial_dimensions) { // ==== 2D ====
   
    int ir,iz,ip;
    int irzp,iirzp;           //AWS20110126

    //Gulli20070810
    propel_arrays_initialize = 0;

    INFO("Generating alpha for 2D");

    alpha1_r.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha1_z.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha1_p.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha2_r.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha2_z.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha2_p.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha3_r.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha3_z.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha3_p.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    
    Nr1_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Nz1_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Np1_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Nr2_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Nz2_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Np2_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Nr3_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Nz3_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Np3_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    
    total_source_function.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Work.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
  
    /* may replace valarrays with pointers for OpenMP performance
    particle_cr_density_arr  .resize(particle.n_rgrid * particle.n_zgrid * particle.n_pgrid); //AWS20110126
    total_source_function_arr.resize(particle.n_rgrid * particle.n_zgrid * particle.n_pgrid); //AWS20110126
    alpha1_r_arr             .resize(particle.n_rgrid * particle.n_zgrid * particle.n_pgrid); //AWS20110126
    alpha1_z_arr             .resize(particle.n_rgrid * particle.n_zgrid * particle.n_pgrid); //AWS20110126
    alpha1_p_arr             .resize(particle.n_rgrid * particle.n_zgrid * particle.n_pgrid); //AWS20110126
    alpha2_r_arr             .resize(particle.n_rgrid * particle.n_zgrid * particle.n_pgrid); //AWS20110126
    alpha2_z_arr             .resize(particle.n_rgrid * particle.n_zgrid * particle.n_pgrid); //AWS20110126
    alpha2_p_arr             .resize(particle.n_rgrid * particle.n_zgrid * particle.n_pgrid); //AWS20110126
    alpha3_r_arr             .resize(particle.n_rgrid * particle.n_zgrid * particle.n_pgrid); //AWS20110126
    alpha3_z_arr             .resize(particle.n_rgrid * particle.n_zgrid * particle.n_pgrid); //AWS20110126
    alpha3_p_arr             .resize(particle.n_rgrid * particle.n_zgrid * particle.n_pgrid); //AWS20110126    
    Work_arr                 .resize(particle.n_rgrid * particle.n_zgrid * particle.n_pgrid); //AWS20110126
    */

    if(galdef.solution_method==21)
    {
    particle_cr_density_arr  =new double[particle.n_rgrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110128
    total_source_function_arr=new double[particle.n_rgrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110128
    alpha1_r_arr             =new double[particle.n_rgrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110128
    alpha1_z_arr             =new double[particle.n_rgrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110128
    alpha1_p_arr             =new double[particle.n_rgrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110126
    alpha2_r_arr             =new double[particle.n_rgrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110126
    alpha2_z_arr             =new double[particle.n_rgrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110126
    alpha2_p_arr             =new double[particle.n_rgrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110126
    alpha3_r_arr             =new double[particle.n_rgrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110126
    alpha3_z_arr             =new double[particle.n_rgrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110126
    alpha3_p_arr             =new double[particle.n_rgrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110126    
    Work_arr                 =new double[particle.n_rgrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110126
    }


    alpha1_r = alpha1_z = alpha1_p = 0.;
    alpha2_r = alpha2_z = alpha2_p = 0.;
    alpha3_r = alpha3_z = alpha3_p = 0.;



      // linear arrays for more efficiency in method 2                             AWS20110126
    if(galdef.solution_method==21)
    {
      irzp=0;                                                                    //AWS20110126
      for (ir = 0; ir < particle.n_rgrid  ; ++ir)                                //AWS20110126
      for (iz = 0; iz < particle.n_zgrid  ; ++iz)                                //AWS20110126          
      for (ip = 0; ip < particle.n_pgrid  ; ++ip)                                //AWS20110126                
      {                                                                          //AWS20110126

         particle_cr_density_arr[irzp] =   particle.cr_density.d2[ir][iz].s[ip]; //AWS20110126

         irzp++;                                                                 //AWS20110126
      }                                                                          //AWS20110126
    }// if solution_method==21
    
    // DIFFUSION term: method according to ApJ 509, 212
    
   

    for (ir = 0; ir < particle.n_rgrid; ++ir) {
	  
      for (iz = 0; iz < particle.n_zgrid; ++iz) {

	for (ip = 0; ip < particle.n_pgrid; ++ip) {
		  
     
   //	  alpha1_z.d2[ir][iz].s[ip] = particle.Dxx.d2[ir][iz].s[ip]*pow(particle.dz,-2.);//AWS20131017
	  alpha1_z.d2[ir][iz].s[ip] = particle.Dzz.d2[ir][iz].s[ip]*pow(particle.dz,-2.);//AWS20131017
	  alpha2_r.d2[ir][iz].s[ip] = particle.Dxx.d2[ir][iz].s[ip]*pow(particle.dr,-2.);
		  

	  if (0 == ir) // use: Dxx/(R[i-1/2] dR)*{Ri(U[i+1]-Ui)-R[i-1](Ui-U[i-1])}; R[i-1]=0
	    alpha3_r.d2[ir][iz].s[ip] = alpha2_r.d2[ir][iz].s[ip] *2.;

	  else { // use: Dxx/(Ri dR)*{R[i+1/2](U[i+1]-Ui)-R[i-1/2](Ui-U[i-1])}
	    
	    alpha1_r.d2[ir][iz].s[ip] = alpha2_r.d2[ir][iz].s[ip]*(1. - particle.dr/2./particle.r[ir]);
	    alpha3_r.d2[ir][iz].s[ip] = alpha2_r.d2[ir][iz].s[ip]*(1. + particle.dr/2./particle.r[ir]);
	
	  }


	}   //ip
    
      }   //iz
    
    }   //ir
      
    alpha2_z = alpha1_z;   
    alpha2_z *= 2.;
    alpha3_z = alpha1_z;                                
      
    alpha2_r *= 2.; // IMOS20061030

    alpha1_r *= factor;   
    alpha1_z *= factor;
    alpha2_r *= factor;   
    alpha2_z *= factor;   
    alpha3_r *= factor;   
    alpha3_z *= factor;

    // CONVECTION term: method according to ApJ 509, 212 with corrections  IMOS20010307
    // convection introduced in v38 using imos code with modifications AWS20010325
    
    if (galdef.convection) {
		     

      for (iz = 0; iz < particle.n_zgrid; ++iz) { // numerator = abs convection velocity in cm s^-1
	/* replaced by v_conv from array AWS20131008
	double az1 = (particle.z[iz] > 0.) ? 
                     (galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz-1]))*1.e5/particle.dz/kpc2cm: 0.; //AWS20110330 removed dt

	double az2 = (galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz]  ))*1.e5/particle.dz/kpc2cm;     //AWS20110330 removed dt
	
	double az3 = (particle.z[iz] < 0.) ? 
                     (galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz+1]))*1.e5/particle.dz/kpc2cm: 0.;//AWS20110330 removed dt
	 */

	// should be exactly the same for v0_conv=0

	double az1 = (particle.z[iz] > 0.) ? fabs(particle.v_conv.d2[0][iz-1].s[0])*1.e5/particle.dz/kpc2cm: 0.; //AWS20131008

	double az2 =                         fabs(particle.v_conv.d2[0][iz  ].s[0])*1.e5/particle.dz/kpc2cm;     //AWS20131008
	
	double az3 = (particle.z[iz] < 0.) ? fabs(particle.v_conv.d2[0][iz+1].s[0])*1.e5/particle.dz/kpc2cm: 0.; //AWS20131008

	if (fabs(particle.z[iz]) < 0.5*particle.dz) 
	  az1 = az2 = az3 = 0.;
	
        az1 *= dt; //AWS20110330
        az2 *= dt; //AWS20110330
        az3 *= dt; //AWS20110330

	for (ip = 1; ip < particle.n_pgrid-1; ++ip) {

       	  /* replaced by dvdz_conv from array AWS20131008
	  const double ap2 = particle.p[ip]  /3.*galdef.dvdz_conv/(particle.p[ip+1]-particle.p[ip])*1.e5/kpc2cm;
	  
	  const double ap3 = particle.p[ip+1]/3.*galdef.dvdz_conv/(particle.p[ip+1]-particle.p[ip])*1.e5/kpc2cm;

	  */
	  	  // keep same scheme (use ir=0) for consistency with original, but later can generalize to convection varying with r


          ir=0;                                                                                                                       //AWS20131008
	  const double ap2 = particle.p[ip]  /3.*particle.dvdz_conv.d2[ir][iz].s[ip] / (particle.p[ip+1]-particle.p[ip])*1.e5/kpc2cm; //AWS20131008
	  
	  const double ap3 = particle.p[ip+1]/3.*particle.dvdz_conv.d2[ir][iz].s[ip]/  (particle.p[ip+1]-particle.p[ip])*1.e5/kpc2cm; //AWS20131008


	  for (ir = 0; ir < particle.n_rgrid; ++ir) {

	    alpha1_z.d2[ir][iz].s[ip] += az1;
	    
	    alpha2_z.d2[ir][iz].s[ip] += az2;   
	    alpha2_p.d2[ir][iz].s[ip] += ap2;
            
	    alpha3_z.d2[ir][iz].s[ip] += az3;   
	    alpha3_p.d2[ir][iz].s[ip] += ap3;
	  }   //ir
            
	}   //ip
	

    
	
        ir=0;                  //AWS20131008
        ip=particle.n_pgrid-1; //AWS20131008
        const double ap2 = particle.p[particle.n_pgrid-1]  /3.*particle.dvdz_conv.d2[ir][iz].s[ip]/(particle.p[particle.n_pgrid-1]-particle.p[particle.n_pgrid-2])*1.e5/kpc2cm; //AWS20131008

	for (ir = 0; ir < particle.n_rgrid; ++ir) {

	  alpha1_z.d2[ir][iz].s[particle.n_pgrid-1] += az1;
	  alpha2_z.d2[ir][iz].s[particle.n_pgrid-1] += az2;
	  alpha3_z.d2[ir][iz].s[particle.n_pgrid-1] += az3;

	  alpha2_p.d2[ir][iz].s[particle.n_pgrid-1] += ap2; //AWS20131008

	}   //ir
        
      }   //iz
      
    }

    // DIFFUSIVE REACCELERATION

    if (galdef.diff_reacc > 0) { 
    
      for (ir = 0; ir < particle.n_rgrid; ++ir) {

	for (iz = 0; iz < particle.n_zgrid; ++iz) {

	  for (ip = 1; ip < particle.n_pgrid-1; ++ip) {

	    // alternative scheme #2, most detailed
	    alpha1_p.d2[ir][iz].s[ip] += 
	      (-(particle.Dpp.d2[ir][iz].s[ip] - particle.Dpp.d2[ir][iz].s[ip-1])
	       /(particle.p[ip] - particle.p[ip-1])
	       + 2.*particle.Dpp.d2[ir][iz].s[ip]
	       /(particle.p[ip+1] - particle.p[ip-1])
	       + 2.*particle.Dpp.d2[ir][iz].s[ip-1]/particle.p[ip-1]) 
	      /(particle.p[ip]-particle.p[ip-1]);

	    alpha2_p.d2[ir][iz].s[ip] +=
	      -(particle.Dpp.d2[ir][iz].s[ip]-particle.Dpp.d2[ir][iz].s[ip-1])
	      /pow(particle.p[ip] - particle.p[ip-1], 2.)
	      + 2.*particle.Dpp.d2[ir][iz].s[ip]
	      /(particle.p[ip+1] - particle.p[ip-1])
	      *(1./(particle.p[ip+1] - particle.p[ip])
		+ 1./(particle.p[ip] - particle.p[ip-1]))
	      + 2.*particle.Dpp.d2[ir][iz].s[ip]
	      /(particle.p[ip] - particle.p[ip-1])
	      /particle.p[ip];
 
	    alpha3_p.d2[ir][iz].s[ip] +=
	      2.*particle.Dpp.d2[ir][iz].s[ip]
	      /(particle.p[ip+1] - particle.p[ip-1])
	      /(particle.p[ip+1] - particle.p[ip]);

	  }   //ip
	
	}   //iz
      
      }   //ir
    
    }   // diffusive reacceleration

    // MOMENTUM LOSSES  IMOS20010307 minor change to make as in ApJ 509, 212
    // AWS20010622 correction to alpha2, consistent with ApJ 509, 212 Table 3 
    // AWS20010622 but note signs of alpha2, alpha3 incorrect in ApJ 509, 212 Table 3
    // AWS20010622 code is correct since dpdt = momentum loss rate is code
    
    if (galdef.momentum_losses) {
     
      for (ir = 0; ir < particle.n_rgrid; ++ir) {

	for (iz = 0; iz < particle.n_zgrid; ++iz) {

	  for (ip = 1; ip < particle.n_pgrid-1; ++ip) {

	    alpha2_p.d2[ir][iz].s[ip] += particle.dpdt.d2[ir][iz].s[ip]/(particle.p[ip+1] - particle.p[ip]);
	    
	    alpha3_p.d2[ir][iz].s[ip] += particle.dpdt.d2[ir][iz].s[ip+1]/(particle.p[ip+1] - particle.p[ip]);
              
	  }   //ip
          
	}   //iz
        
      }   //ir 
      
    }   //momentum losses

    // SET VALUES AT EXTREMES OF P GRID
    
    
    for (ir = 0; ir < particle.n_rgrid; ++ir) {

      for (iz = 0; iz < particle.n_zgrid; ++iz) {

	alpha1_p.d2[ir][iz].s[0] = alpha1_p.d2[ir][iz].s[1];
	alpha2_p.d2[ir][iz].s[0] = alpha2_p.d2[ir][iz].s[1];
	alpha3_p.d2[ir][iz].s[0] = alpha3_p.d2[ir][iz].s[1];
	
	alpha1_p.d2[ir][iz].s[particle.n_pgrid-1] = alpha1_p.d2[ir][iz].s[particle.n_pgrid-2];
	alpha2_p.d2[ir][iz].s[particle.n_pgrid-1] = alpha2_p.d2[ir][iz].s[particle.n_pgrid-2];
	alpha3_p.d2[ir][iz].s[particle.n_pgrid-1] = alpha3_p.d2[ir][iz].s[particle.n_pgrid-2];
    
      }//iz
    
    }//ir 

    alpha1_p *= dt;
    alpha2_p *= dt;
    alpha3_p *= dt;
    
    double f_use = galdef.prop_r + galdef.prop_z + galdef.prop_p; // used only by method = 1
    

    /* change to method-independent alpha definitions without f_use AWS20110220
    // FRAGMENTATION
    
    if (galdef.fragmentation) {
       
      Work = particle.fragment;         // use Work to avoid memory leak
      Work*= (dt/f_use);
      alpha2_r += Work;                  // fragment used f_use times
      alpha2_z += Work;                         
      alpha2_p += Work;                          
 
    }

    // DECAY
    
    if (particle.t_half != 0 && galdef.radioactive_decay) {

      ostringstream rBuf;
      rBuf << "Radioactive decay of " << particle.name << " with half life " << particle.t_half;
      INFO(rBuf.str());
 
      Work = particle.decay;            // use Work to avoid memory leak
      Work *= (dt/f_use);
      alpha2_r += Work;                  // decay used f_use times
      alpha2_z += Work;                         
      alpha2_p += Work;
 
    }

    */



    //AWS20110120
    //method-independent alpha definitions without f_use, using alpha_2 only once AWS20110220
    // assign to exactly one of the dimensions, preferring momentum if present, otherwise R or z
    // FRAGMENTATION
    
    if (galdef.fragmentation) {
       
      Work  = particle.fragment;         // use Work to avoid memory leak
      Work *= dt;

      if      (galdef.prop_p)alpha2_p += Work; 
      else if (galdef.prop_r)alpha2_r += Work;                  
      else if (galdef.prop_z)alpha2_z += Work;                         
                          
 
    }

    // DECAY
    
    if (particle.t_half != 0 && galdef.radioactive_decay) {

      ostringstream rBuf;
      rBuf << "Radioactive decay of " << particle.name << " with half life " << particle.t_half;
      INFO(rBuf.str());
 
      Work  = particle.decay;            // use Work to avoid memory leak
      Work *= dt;

      if      (galdef.prop_p)alpha2_p += Work; 
      else if (galdef.prop_r)alpha2_r += Work;                  
      else if (galdef.prop_z)alpha2_z += Work; 
 
    }

    buf.str("");
    buf<<"alpha1_r:"<< alpha1_r.d2[0][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha1_r.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha2_r.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha3_r.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha1_z.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha2_z.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha3_z.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha1_p.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha2_p.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha3_p.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    
    if (galdef.verbose==10) {
       
      cout<<"alpha1_r:"<<endl;  alpha1_r.print();
      cout<<"alpha2_r:"<<endl;  alpha2_r.print();
      cout<<"alpha3_r:"<<endl;  alpha3_r.print();
      cout<<"alpha1_z:"<<endl;  alpha1_z.print();
      cout<<"alpha2_z:"<<endl;  alpha2_z.print();
      cout<<"alpha3_z:"<<endl;  alpha3_z.print();
      cout<<"alpha1_p:"<<endl;  alpha1_p.print();
      cout<<"alpha2_p:"<<endl;  alpha2_p.print();
      cout<<"alpha3_p:"<<endl;  alpha3_p.print();
    
    }

    ostringstream tBuf;
    tBuf << "galdef_ID " << galdef.galdef_ID << " start timestep dt = " << dt << " " << particle.name << endl;
    DEBUGLOG(tBuf.str());
    
    int timestep_mode;
    int timestep_mode2done = 0;
    
    for (timestep_mode = 1; timestep_mode <= 2; ++timestep_mode) {
      
      buf.str("");
      buf<<"    timestep_mode="<<timestep_mode;
      INFO(buf.str());
      t = 0;
         
      while ((1 == start) || 
	     (dt > end_timestep_sec/galdef.timestep_factor*0.9999) || 
	     (2 == timestep_mode && 0 == timestep_mode2done)) {
         
	if (0 == start && 1 == timestep_mode) {


	  //    while (1 == start || dt > end_timestep_sec/galdef.timestep_factor*0.9999) {   //AWS20110119
      
	  //      if (0 == start) {                                                           //AWS20110119

	dt *= galdef.timestep_factor;
	tBuf.str("");
	tBuf<<" galdef_ID "<<galdef.galdef_ID<<" new timestep dt="<<dt<<"  "
	   <<particle.name;
	DEBUGLOG(tBuf.str());

	alpha1_r *= galdef.timestep_factor;
	alpha2_r *= galdef.timestep_factor;
	alpha3_r *= galdef.timestep_factor;
	alpha1_z *= galdef.timestep_factor;
	alpha2_z *= galdef.timestep_factor;
	alpha3_z *= galdef.timestep_factor;
	alpha1_p *= galdef.timestep_factor;
	alpha2_p *= galdef.timestep_factor;
	alpha3_p *= galdef.timestep_factor;
      
      }
      
      start = 0;
      
      Nr1_ = alpha1_r;   Nr1_*= -0.5;
      Nr2_ = alpha2_r;   Nr2_*=  0.5;   Nr2_ += 1.0;   
      Nr3_ = alpha3_r;   Nr3_*= -0.5;
      Nz1_ = alpha1_z;   Nz1_*= -0.5;
      Nz2_ = alpha2_z;   Nz2_*=  0.5;   Nz2_ += 1.0;   
      Nz3_ = alpha3_z;   Nz3_*= -0.5;
      Np1_ = alpha1_p;   Np1_*= -0.5;
      Np2_ = alpha2_p;   Np2_*=  0.5;   Np2_ += 1.0;   
      Np3_ = alpha3_p;   Np3_*= -0.5;
      
      buf.str("");
      buf<<"Nr1_    ;"<< Nr1_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Nr2_    :"<< Nr2_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Nr3_    :"<< Nr3_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Nz1_    :"<< Nz1_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Nz2_    :"<< Nz2_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Nz3_    :"<< Nz3_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Np1_    :"<< Np1_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Np2_    :"<< Np2_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Np3_    :"<< Np3_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      
      total_source_function = particle.primary_source_function + particle.secondary_source_function;
      //total_source_function += particle.secondary_source_function;

      total_source_function *= dt; //AWS20110131 put dt here for consistency with other terms

      if(galdef.solution_method==21)
      {
      // valarrays for more efficiency in method 2                                 AWS20110126
      irzp=0;                                                                    //AWS20110126
      for (ir = 0; ir < particle.n_rgrid  ; ++ir)                                //AWS20110126
      for (iz = 0; iz < particle.n_zgrid  ; ++iz)                                //AWS20110126          
      for (ip = 0; ip < particle.n_pgrid  ; ++ip)                                //AWS20110126                
      {                                                                          //AWS20110126

       total_source_function_arr[irzp] = total_source_function.d2[ir][iz].s[ip]; //AWS20110126

                    alpha1_r_arr[irzp] =              alpha1_r.d2[ir][iz].s[ip]; //AWS20110126
                    alpha1_z_arr[irzp] =              alpha1_z.d2[ir][iz].s[ip]; //AWS20110126
                    alpha1_p_arr[irzp] =              alpha1_p.d2[ir][iz].s[ip]; //AWS20110126
                    alpha2_r_arr[irzp] =              alpha2_r.d2[ir][iz].s[ip]; //AWS20110126
                    alpha2_z_arr[irzp] =              alpha2_z.d2[ir][iz].s[ip]; //AWS20110126
                    alpha2_p_arr[irzp] =              alpha2_p.d2[ir][iz].s[ip]; //AWS20110126
                    alpha3_r_arr[irzp] =              alpha3_r.d2[ir][iz].s[ip]; //AWS20110126
                    alpha3_z_arr[irzp] =              alpha3_z.d2[ir][iz].s[ip]; //AWS20110126
                    alpha3_p_arr[irzp] =              alpha3_p.d2[ir][iz].s[ip]; //AWS20110126

       irzp++;                                                                   //AWS20110126
      }                                                                          //AWS20110126
      }// if solution_method==21

      if (0 == total_source_function.max()) 
	return 0;

      // PROPAGATION iterations
      
      //     int  irept;                                                      //AWS20110119
      //     for (irept = 1; irept <= galdef.timestep_repeat; ++irept) {      //AWS20110119
	
	int nrept;                                                            //AWS20111019
	if (1 == timestep_mode)                                               //AWS20111019
	  nrept = galdef.timestep_repeat;                                     //AWS20111019

	if (2 == timestep_mode)                                               //AWS20111019
	  nrept = galdef.timestep_repeat2; //AWS20111019
	
	if (1 == timestep_mode)
        {
          method = 1;  // fixed to ensure stability
	  buf.str(""); buf<<"timestep_mode 1: using solution method "<<method<< " timestep="<<dt/year2sec<< " yr";
          INFO(buf);       
        }

	if (2 == timestep_mode)
        {
          if(galdef.solution_method== 1){method=1;         }    //AWS20110131
	  if(galdef.solution_method== 2){method=2; turbo=0;}    //AWS20110131
          if(galdef.solution_method==21){method=2; turbo=1;}    //AWS20110131

          if(method==1){buf.str(""); buf<<"timestep_mode 2: using  method "<<method<<" Crank-Nicolson"
                                                          << " timestep="<<dt/year2sec<< " yr";}
	  if(method==2){buf.str(""); buf<<"timestep_mode 2: using  method "<<method<<" fully time explicit, "
		                        <<" turbo="<<turbo<< " timestep="<<dt/year2sec<< " yr" ;}
          INFO(buf);       

	  // following required to transfer from method=1 in timestep_mode 1
         if(galdef.solution_method==21)
         {
          irzp=0;                                                                    //AWS20110130
          for (ir = 0; ir < particle.n_rgrid  ; ++ir)                                //AWS20110130
          for (iz = 0; iz < particle.n_zgrid  ; ++iz)                                //AWS20110130          
          for (ip = 0; ip < particle.n_pgrid  ; ++ip)                                //AWS20110130                
          {                                                                          //AWS20110130
           particle_cr_density_arr[irzp] =   particle.cr_density.d2[ir][iz].s[ip];   //AWS20110130
           irzp++;                                                                   //AWS20110130
          }                                                                          //AWS20110130
         }// if solution_method==21


        }// timestep_mode==2

        int  irept;                                                           //AWS20111018
	for (irept = 1; irept <= nrept; ++irept) {                            //AWS20111019

	if (method == 1) // Crank-Nicolson method AWS20111018
	{


	// Z propagation
	//Gulli20070821
#pragma omp parallel default(shared) private(iz,ir,ip)
	if (1 == galdef.prop_z) {



	  valarray<double> Nz0(0., particle.n_zgrid);
	  valarray<double> Nz1(0., particle.n_zgrid);
	  valarray<double> Nz2(0., particle.n_zgrid);
	  valarray<double> Nz3(0., particle.n_zgrid);
	  valarray<double> Rz (0., particle.n_zgrid);

	  //Gulli20070821
#pragma omp for schedule(runtime) 
	  for (ir = 0; ir < particle.n_rgrid; ++ir) {
	    
	    for (ip = 0; ip < particle.n_pgrid; ++ip) {
	      
	      for (iz = 0; iz < particle.n_zgrid; ++iz) {
                  
		Nz1[iz] = Nz1_.d2[ir][iz].s[ip];
		Nz2[iz] = Nz2_.d2[ir][iz].s[ip];
		Nz3[iz] = Nz3_.d2[ir][iz].s[ip];
		
		//		Nz0[iz] = total_source_function.d2[ir][iz].s[ip]*dt/f_use;
		Nz0[iz] = total_source_function.d2[ir][iz].s[ip]  /f_use;//AWS20110131
		Nz0[iz] += (1. - 0.5*alpha2_z.d2[ir][iz].s[ip])*particle.cr_density.d2[ir][iz].s[ip];
             
		if (iz > 0) 
		  Nz0[iz] += 0.5*alpha1_z.d2[ir][iz].s[ip] 
		    *particle.cr_density.d2[ir][iz-1].s[ip];
                    
		if (iz < particle.n_zgrid-1) 
		  Nz0[iz] += 0.5*alpha3_z.d2[ir][iz].s[ip] 
		    *particle.cr_density.d2[ir][iz+1].s[ip];
                    
	      } 

	      tridag(&Nz1[0], &Nz2[0], &Nz3[0], &Nz0[0], &Rz[0], particle.n_zgrid);  

	      for (iz = 0; iz < particle.n_zgrid; ++iz)
		particle.cr_density.d2[ir][iz].s[ip] = (Rz[iz] < 0.) ? 0. : Rz[iz]; //IMOS20030217

	     if (galdef.spatial_bound_conds==1) //AWS20131028
	     {
	      iz = 0;                                     // boundary condition
	      particle.cr_density.d2[ir][iz].s[ip] = 0;

	      iz = particle.n_zgrid-1;                    // boundary condition
	      particle.cr_density.d2[ir][iz].s[ip] = 0;
             }
     
	    }   //ip
            
	  }   //ir
      

	
	}   //prop_z

	// P propagation

	//Gulli20070821
#pragma omp parallel default(shared) private(iz,ir,ip)
	if (1 == galdef.prop_p) { 
	  

	  
	  valarray<double> Np0(0., particle.n_pgrid);
	  valarray<double> Np1(0., particle.n_pgrid);
	  valarray<double> Np2(0., particle.n_pgrid);
	  valarray<double> Np3(0., particle.n_pgrid);
	  valarray<double> Rp (0., particle.n_pgrid);

	  //Gulli20070821
#pragma omp for schedule(runtime) 
	  for (ir = 0; ir < particle.n_rgrid; ++ir) {
	    
	    for (iz = 0; iz < particle.n_zgrid; ++iz) {
	      
	      for (ip = 0; ip < particle.n_pgrid; ++ip) {

		Np1[ip] = Np1_.d2[ir][iz].s[ip];
		Np2[ip] = Np2_.d2[ir][iz].s[ip];
		Np3[ip] = Np3_.d2[ir][iz].s[ip];
		
		//		Np0[ip] = total_source_function.d2[ir][iz].s[ip]*dt/f_use;
		Np0[ip] = total_source_function.d2[ir][iz].s[ip]   /f_use;//AWS20110131
                        
		Np0[ip] += (1. - 0.5*alpha2_p.d2[ir][iz].s[ip])
		  *particle.cr_density.d2[ir] [iz].s[ip];
                        
		if (ip > 0) 
		  Np0[ip] += 0.5*alpha1_p.d2[ir][iz].s[ip] 
		    *particle.cr_density.d2[ir][iz].s[ip-1];
                    
		if (ip < particle.n_pgrid-1) 
		  Np0[ip] += 0.5*alpha3_p.d2[ir][iz].s[ip] 
		    *particle.cr_density.d2[ir][iz].s[ip+1];
                     
	      } 

	      tridag(&Np1[0], &Np2[0], &Np3[0], &Np0[0], &Rp[0], particle.n_pgrid);  

	      for (ip = 0; ip < particle.n_pgrid; ++ip)
		particle.cr_density.d2[ir][iz].s[ip] = (Rp[ip] < 0.) ? 0.: Rp[ip]; //IMOS20030217
	      
	      // ir=particle.n_rgrid-1;  particle.cr_density.d2[ir][iz].s[ip]=0.0; // boundary condition
	      // if(irept%  1==0) cout<<"\nparticle.cr_density after iteration "<<irept<<endl;
	      // particle.cr_density.print();
	    }   //iz
               
	  }   //ir
	  

	
	}   //prop_p

	// R propagation
	
	//Gulli20070821
#pragma omp parallel default(shared) private(iz,ir,ip)
	if (1 == galdef.prop_r) {



	  valarray<double> Nr0(0., particle.n_rgrid);
	  valarray<double> Nr1(0., particle.n_rgrid);
	  valarray<double> Nr2(0., particle.n_rgrid);
	  valarray<double> Nr3(0., particle.n_rgrid);
	  valarray<double> Rr (0., particle.n_rgrid);

 //Gulli20070821
#pragma omp for schedule(runtime) 
	  for (iz = 0; iz < particle.n_zgrid; ++iz) {
	    
	    for (ip = 0; ip < particle.n_pgrid; ++ip) {

	      for (ir = 0; ir < particle.n_rgrid; ++ir) {

		Nr1[ir] = Nr1_.d2[ir][iz].s[ip];
		Nr2[ir] = Nr2_.d2[ir][iz].s[ip];
		Nr3[ir] = Nr3_.d2[ir][iz].s[ip];
		
		//		Nr0[ir] = total_source_function.d2[ir][iz].s[ip]*dt/f_use;
		Nr0[ir] = total_source_function.d2[ir][iz].s[ip]   /f_use; //AWS20110131
	
		Nr0[ir] += (1. - 0.5*alpha2_r.d2[ir][iz].s[ip])
		  *particle.cr_density.d2[ir][iz].s[ip];
                        
		if (ir > 0) 
		  Nr0[ir] += 0.5*alpha1_r.d2[ir][iz].s[ip] 
		    *particle.cr_density.d2[ir-1][iz].s[ip];
                    
		if (ir < particle.n_rgrid-1) 
		  Nr0[ir] += 0.5*alpha3_r.d2[ir][iz].s[ip] 
		    *particle.cr_density.d2[ir+1][iz].s[ip];
                    
	      }

	      tridag(&Nr1[0], &Nr2[0], &Nr3[0], &Nr0[0], &Rr[0], particle.n_rgrid);  

	      for (ir = 0; ir < particle.n_rgrid; ++ir)
		particle.cr_density.d2[ir][iz].s[ip] = (Rr[ir] < 0.) ? 0.: Rr[ir]; //IMOS20030217
  
	     if( galdef.spatial_bound_conds==1) //AWS20131028
	     { 
	      ir = particle.n_rgrid-1;                  // boundary condition
	      particle.cr_density.d2[ir][iz].s[ip] = 0;
	     }

	    }   //ip
	  
	  }   //iz
                 
	}   //prop_r



	  // apply spatial boundary conditions AWS20110130
          // since the boundaries are not completely set in the scheme above
          
	 if( galdef.spatial_bound_conds==1) //AWS20131028
	 { 
          // in R, only at outer boundary
	  // NB ir must not be private or wrong result

	   ir = particle.n_rgrid-1;
#pragma omp parallel for schedule(runtime)  private(iz,ip) default(shared)	    
	   for (iz = 0; iz < particle.n_zgrid; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid; ++ip)
            {
	       particle.cr_density.d2[ir][iz].s[ip] = 0.;
	    }//ip
	   }//iz

	   // in z at both boundaries
	  // NB iz must not be private or wrong result

	   iz = 0;
#pragma omp parallel for schedule(runtime)  private(ir,ip) default(shared)	    
	   for (ir = 0; ir < particle.n_rgrid; ++ir)
           {
	    for (ip = 0; ip < particle.n_pgrid; ++ip)
            {
	       particle.cr_density.d2[ir][iz].s[ip] = 0.;
	    }//ip
	   }//ir

	   iz = particle.n_zgrid-1;
#pragma omp parallel for schedule(runtime)  private(ir,ip) default(shared)	    
	   for (ir = 0; ir < particle.n_rgrid; ++ir)
           {
	    for (ip = 0; ip < particle.n_pgrid; ++ip)
            {
	       particle.cr_density.d2[ir][iz].s[ip] = 0.;
	    }//ip
	   }//ir

	 } // if spatial_bound_conds==1

	// end of updating step for Crank-Nicolson method
       } //method == 1 AWS20110118



	//--------------------------------------------------


	if (method == 2) // explicit updating method AWS20110118
	{
         
	  // group by boundaries for efficient optimization avoiding conditionals
	  // loop indices are private by default so no need to specify, but perhaps need to for inner loops, so do it.



	  int nrzp =  particle.n_rgrid *  particle.n_zgrid *  particle.n_pgrid;

          if (turbo==0)
	  { 
#pragma omp parallel for schedule(runtime)  private(ir,iz,ip) default(shared)
	  for (ir = 0; ir < particle.n_rgrid  ; ++ir)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {

	      Work.d2[ir][iz].s[ip]  =                           total_source_function.d2[ir  ][iz  ].s[ip];//AWS20110131

	      Work.d2[ir][iz].s[ip] -= alpha2_r.d2[ir][iz].s[ip] * particle.cr_density.d2[ir  ][iz  ].s[ip];
	      Work.d2[ir][iz].s[ip] -= alpha2_z.d2[ir][iz].s[ip] * particle.cr_density.d2[ir  ][iz  ].s[ip];  
	      Work.d2[ir][iz].s[ip] -= alpha2_p.d2[ir][iz].s[ip] * particle.cr_density.d2[ir  ][iz  ].s[ip];
	      

	      /* to removed to assist optimization/parallelization
              if(0==1)
	      {
              cout<<"method 2: dt="<<dt;
              cout<<"          total_source_function*dt="<<total_source_function.d2[ir][iz].s[ip]*dt<<" " ;
              cout<<"          alpha1_r="<<alpha1_r.d2[ir][iz].s[ip]<<" " ;
              cout<<"          alpha2_r="<<alpha2_r.d2[ir][iz].s[ip]<<" " ;
              cout<<"          alpha3_r="<<alpha3_r.d2[ir][iz].s[ip]<<" " ;
              cout<<"particle.cr_density="<<particle.cr_density.d2[ir][iz].s[ip]<<" " ;
              cout<<"          Work=    "<<    Work.d2[ir][iz].s[ip]<<endl;
	      }
	      */

	    }//ip
	   }//iz
	  }//ir
	 


	  
	 
#pragma omp parallel for schedule(runtime)  private(ir,iz,ip) default(shared)
	  for (ir = 1; ir < particle.n_rgrid  ; ++ir)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
	      Work.d2[ir][iz].s[ip] += alpha1_r.d2[ir][iz].s[ip] * particle.cr_density.d2[ir-1][iz  ].s[ip];
	    }//ip
	   }//iz
	  }//ir
	  






       
       
#pragma omp parallel for schedule(runtime)  private(ir,iz,ip) default(shared)
	  for (ir = 0; ir < particle.n_rgrid-1; ++ir)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
	      Work.d2[ir][iz].s[ip] += alpha3_r.d2[ir][iz].s[ip] * particle.cr_density.d2[ir+1][iz  ].s[ip];
	    }//ip
	   }//iz
	  }//ir
       



       
	

#pragma omp parallel for schedule(runtime)  private(ir,iz,ip) default(shared)
	  for (ir = 0; ir < particle.n_rgrid  ; ++ir)
          {
	   for (iz = 1; iz < particle.n_zgrid  ; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {	                  
	      Work.d2[ir][iz].s[ip] += alpha1_z.d2[ir][iz].s[ip] * particle.cr_density.d2[ir  ][iz-1].s[ip];
	    }//ip
	   }//iz
	  }//ir


#pragma omp parallel for schedule(runtime)  private(ir,iz,ip) default(shared)
	  for (ir = 0; ir < particle.n_rgrid  ; ++ir)
          {
	   for (iz = 0; iz < particle.n_zgrid-1; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {	                  
	      Work.d2[ir][iz].s[ip] += alpha3_z.d2[ir][iz].s[ip] * particle.cr_density.d2[ir  ][iz+1].s[ip];
	    }//ip
	   }//iz
	  }//ir

	 


	
	 

#pragma omp parallel for schedule(runtime)  private(ir,iz,ip) default(shared)
	  for (ir = 0; ir < particle.n_rgrid  ; ++ir)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	    for (ip = 1; ip < particle.n_pgrid  ; ++ip)
            {
	      Work.d2[ir][iz].s[ip] += alpha1_p.d2[ir][iz].s[ip] * particle.cr_density.d2[ir  ][iz  ].s[ip-1];
	    }//ip
	   }//iz
	  }//ir

#pragma omp parallel for schedule(runtime)  private(ir,iz,ip) default(shared)
	  for (ir = 0; ir < particle.n_rgrid  ; ++ir)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid-1; ++ip)
            {
	      Work.d2[ir][iz].s[ip] += alpha3_p.d2[ir][iz].s[ip] * particle.cr_density.d2[ir  ][iz  ].s[ip+1];
	    }//ip
	   }//iz
	  }//ir

	 



	 
       
#pragma omp parallel for schedule(runtime)  private(ir,iz,ip) default(shared)
	  for (ir = 0; ir < particle.n_rgrid  ; ++ir) //AWS20110126
          {
           for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
	      particle.cr_density.d2[ir][iz].s[ip] += Work.d2[ir][iz].s[ip] ;
	    }//ip
	   }//iz
	  }//ir



	  // apply spatial boundary conditions
          if( galdef.spatial_bound_conds==1) //AWS20131028
	  { 

          // in R, only at outer boundary
	  // NB ir must not be private or wrong result

	   ir = particle.n_rgrid-1;
#pragma omp parallel for schedule(runtime)  private(iz,ip) default(shared)	    
	   for (iz = 0; iz < particle.n_zgrid; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid; ++ip)
            {
	       particle.cr_density.d2[ir][iz].s[ip] = 0.;
	    }//ip
	   }//iz

	   // in z at both boundaries
	  // NB iz must not be private or wrong result

	   iz = 0;
#pragma omp parallel for schedule(runtime)  private(ir,ip) default(shared)	    
	   for (ir = 0; ir < particle.n_rgrid; ++ir)
           {
	    for (ip = 0; ip < particle.n_pgrid; ++ip)
            {
	       particle.cr_density.d2[ir][iz].s[ip] = 0.;
	    }//ip
	   }//ir

	   iz = particle.n_zgrid-1;
#pragma omp parallel for schedule(runtime)  private(ir,ip) default(shared)	    
	   for (ir = 0; ir < particle.n_rgrid; ++ir)
           {
	    for (ip = 0; ip < particle.n_pgrid; ++ip)
            {
	       particle.cr_density.d2[ir][iz].s[ip] = 0.;
	    }//ip
	   }//ir
	
	  } // if spatial_bound_conds==1

	 }//turbo==0


	  //----------------------------------------------------------------------------


	  if (turbo==1) //AWS20110126 
	 { 
	  
	   #pragma omp parallel for schedule(runtime)   default(shared)
	    for (irzp=0; irzp<nrzp  ; ++irzp)
	    {
              /*
	      Work_arr[irzp]         =                           total_source_function_arr[irzp]; //AWS20110131
	      
	      
	      Work_arr[irzp]        -= alpha2_r_arr[irzp]        * particle_cr_density_arr[irzp];             
	      Work_arr[irzp]        -= alpha2_z_arr[irzp]        * particle_cr_density_arr[irzp];             
	      Work_arr[irzp]        -= alpha2_p_arr[irzp]        * particle_cr_density_arr[irzp];     
	      */

	      Work_arr[irzp]        =  total_source_function_arr[irzp]                                            //AWS20110201
                                                  -(alpha2_r_arr[irzp]+alpha2_z_arr[irzp]+ alpha2_p_arr[irzp])    //AWS20110201
                                                  *                             particle_cr_density_arr[irzp];    //AWS20110201
   
      
	    }//irzp




	   

	  // r direction -----------------------------------------------------------


	  
	  
	    #pragma omp parallel for schedule(runtime)  private(ir,iz,ip,irzp,iirzp) default(shared) 
	   
	  for (ir = 1; ir < particle.n_rgrid  ; ++ir)
          {
  	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
             irzp =  particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid * iz;                       //AWS20110201
	    iirzp =  particle.n_zgrid * particle.n_pgrid * (ir-1) +  particle.n_pgrid * iz;                       //AWS20110201
	    
            for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
	      //               irzp =  particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid * iz + ip;//AWS20110201
	      //	      iirzp =  particle.n_zgrid * particle.n_pgrid * (ir-1) +  particle.n_pgrid * iz + ip;//AWS20110201
     
	      Work_arr[irzp]        += alpha1_r_arr[irzp]        * particle_cr_density_arr[iirzp];

               irzp++;                                                                                            //AWS20110201
	      iirzp++;                                                                                            //AWS20110201

	    }//ip               	    
	   }//iz
 	  }//ir



	  #pragma omp parallel for schedule(runtime)  private(ir,iz,ip,irzp,iirzp) default(shared) 
      
	  for (ir = 0; ir < particle.n_rgrid-1; ++ir)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
             irzp =  particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid * iz;                       //AWS20110201
	    iirzp =  particle.n_zgrid * particle.n_pgrid * (ir+1) +  particle.n_pgrid * iz;                       //AWS20110201

	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
	      //               irzp =  particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid * iz + ip;
	      //	      iirzp =  particle.n_zgrid * particle.n_pgrid * (ir+1) +  particle.n_pgrid * iz + ip;
     
	      Work_arr[irzp]        += alpha3_r_arr[irzp]        * particle_cr_density_arr[iirzp];


               irzp++;                                                                                            //AWS20110201
	      iirzp++;                                                                                            //AWS20110201

	    }//ip
	   }//iz
	  }//ir



	  


	  // z direction ---------------------------------------------------------------





       
	  #pragma omp parallel for schedule(runtime)  private(ir,iz,ip,irzp,iirzp) default(shared) 
	  for (ir = 0; ir < particle.n_rgrid  ; ++ir)
          {
	   for (iz = 1; iz < particle.n_zgrid  ; ++iz)
           {
             irzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid *  iz        ;               //AWS20110201
	    iirzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid * (iz-1)     ;               //AWS20110201

	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
	      //                irzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid *  iz    + ip;
	      //	       iirzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid * (iz-1) + ip;
     
	      Work_arr[irzp]        += alpha1_z_arr[irzp]        * particle_cr_density_arr[iirzp];


               irzp++;                                                                                            //AWS20110201
	      iirzp++;                                                                                            //AWS20110201

	    }//ip
	   }//iz
	  }//ir



	  #pragma omp parallel for schedule(runtime)  private(ir,iz,ip,irzp,iirzp) default(shared) 
	  for (ir = 0; ir < particle.n_rgrid  ; ++ir)
          {
	   for (iz = 0; iz < particle.n_zgrid-1; ++iz)
           {
             irzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid *  iz        ;               //AWS20110201
	    iirzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid * (iz+1)     ;               //AWS20110201

	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
	      //             irzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid *  iz    + ip;
	      //	    iirzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid * (iz+1) + ip;
     
	      Work_arr[irzp]        += alpha3_z_arr[irzp]        * particle_cr_density_arr[iirzp];

               irzp++;                                                                                            //AWS20110201
	      iirzp++;                                                                                            //AWS20110201

	    }//ip
	   }//iz
	  }//ir



	 


	  // p direction ---------------------------------------------------------------




       
	
	  #pragma omp parallel for schedule(runtime)  private(ir,iz,ip,irzp,iirzp) default(shared) 
	  for (ir = 0; ir < particle.n_rgrid  ; ++ir)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
             irzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid *  iz     ;                  //AWS20110201
	    iirzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid *  iz  - 1;                  //AWS20110201

	    for (ip = 1; ip < particle.n_pgrid  ; ++ip)
            {
	      //            irzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid *  iz    + ip;
	      //	   iirzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid *  iz    + ip - 1;
     
	      Work_arr[irzp]        += alpha1_p_arr[irzp]        * particle_cr_density_arr[iirzp];

               irzp++;                                                                                            //AWS20110201
	      iirzp++;                                                                                            //AWS20110201

	    }//ip
	   }//iz
	  }//ir



	  #pragma omp parallel for schedule(runtime)  private(ir,iz,ip,irzp,iirzp) default(shared) 
	  for (ir = 0; ir < particle.n_rgrid  ; ++ir)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
             irzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid *  iz     ;                  //AWS20110201
	    iirzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid *  iz  + 1;                  //AWS20110201

	    for (ip = 0; ip < particle.n_pgrid-1; ++ip)
            {
	      //             irzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid *  iz    + ip;
	      //	    iirzp = particle.n_zgrid * particle.n_pgrid *  ir    +  particle.n_pgrid *  iz    + ip + 1;
     
	      Work_arr[irzp]        += alpha3_p_arr[irzp]        * particle_cr_density_arr[iirzp];

               irzp++;                                                                                            //AWS20110201
	      iirzp++;                                                                                            //AWS20110201

	    }//ip
	   }//iz
	  }//ir



	   #pragma omp parallel for schedule(runtime)   default(shared)
	    for (irzp=0; irzp<nrzp  ; ++irzp)                                                                     //AWS20110201
	    {
	     particle_cr_density_arr[irzp]        +=     Work_arr[irzp];
	    }


	    /* AWS20110201
	  #pragma omp parallel for schedule(runtime)   default(shared) private(iz,ip,irzp)	  
	  for (ir = 0; ir < particle.n_rgrid  ; ++ir)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
              irzp =  particle.n_zgrid * particle.n_pgrid * ir +  particle.n_pgrid * iz + ip;
	   
	      particle_cr_density_arr[irzp]        +=     Work_arr[irzp];

	    }//ip
	   }//iz
	  }//ir
	    */
	 	 

	


	  // apply spatial boundary conditions
          
	 if( galdef.spatial_bound_conds==1) //AWS20131028
	 { 
          // in R, only at outer boundary
	  // NB ir must not be private or wrong result

	   ir = particle.n_rgrid-1;
#pragma omp parallel for schedule(runtime)  private(iz,ip,irzp) default(shared)	    
	   for (iz = 0; iz < particle.n_zgrid; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid; ++ip)
            {
             irzp =  particle.n_zgrid * particle.n_pgrid * ir +  particle.n_pgrid * iz + ip;
	     particle_cr_density_arr[irzp]        = 0.;
	    }//ip
	   }//iz

	   // in z at both boundaries
	  // NB iz must not be private or wrong result

	   iz = 0;
#pragma omp parallel for schedule(runtime)  private(ir,ip,irzp) default(shared)	    
	   for (ir = 0; ir < particle.n_rgrid; ++ir)
           {
	    for (ip = 0; ip < particle.n_pgrid; ++ip)
            {
             irzp =  particle.n_zgrid * particle.n_pgrid * ir +  particle.n_pgrid * iz + ip;
	     particle_cr_density_arr[irzp] = 0.;
	    }//ip
	   }//ir

	   iz = particle.n_zgrid-1;
#pragma omp parallel for schedule(runtime)  private(ir,ip,irzp) default(shared)	    
	   for (ir = 0; ir < particle.n_rgrid; ++ir)
           {
	    for (ip = 0; ip < particle.n_pgrid; ++ip)
            {
             irzp =  particle.n_zgrid * particle.n_pgrid * ir +  particle.n_pgrid * iz + ip;
	     particle_cr_density_arr[irzp] = 0.;
	    }//ip
	   }//ir

	 } //if spatial_bound_conds==1	

	   
// following needed only for diagnostics, so conditional for efficiency. 

	if(irept%galdef.timestep_diagnostics==0)
#pragma omp parallel for schedule(runtime)   default(shared) private(ir,iz,ip,irzp)	  
	  for (ir = 0; ir < particle.n_rgrid  ; ++ir)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
              irzp =  particle.n_zgrid * particle.n_pgrid * ir +  particle.n_pgrid * iz + ip;
	   
	      particle.cr_density.d2[ir][iz].s[ip] =      particle_cr_density_arr[irzp];                    
	      
	    
	    }//ip
	   }//iz
	  }//ir

	 }//turbo==1


	}// method == 2 



	//------------------------------------------------------------------


	if (0 == irept%galdef.timestep_print) {

	  ostringstream repBuf;
	  repBuf << " particle.cr_density after iteration "<<irept<<" with time step " <<dt/year2sec<<" yr";
	  INFO(repBuf.str());
	  particle.cr_density.print();   cout<<endl;
            
	}
        

	buf.str("");
	buf<<"maximum="<<particle.cr_density.max();
	buf<<"    propel iteration "<<irept<<endl;     

	if(irept%galdef.timestep_diagnostics==0)
	{
          converged =                                                //AWS20110114
	  propel_diagnostics (particle, alpha1_r, alpha1_z, alpha1_p, 
			      alpha2_r, alpha2_z, alpha2_p,
			      alpha3_r, alpha3_z, alpha3_p,
			      total_source_function, dt);


         if(galdef.solution_convergence == 1 && converged == 1) //AWS20111020
	 {
	  buf.str("");
          buf<<" solution converged at irept="<<irept<<" exiting this timestep";
          INFO(buf.str());
	 
          break;// out of irept loop
         }
        }// if irept 
         

      }   //irept
    
	if (2 == timestep_mode)                    //AWS20111019
	  timestep_mode2done = 1;                  //AWS20111019
      
      }  //  dt>=end_timestep_sec
    
    }  //  timestep_mode                           //AWS20111019

   
    // do this here since method 2 turbo does not do it in loop (for efficiency) except when diagnostics required

    if (galdef.solution_method == 21 && galdef.timestep_repeat2 > 0) //AWS20110130
	  for (ir = 0; ir < particle.n_rgrid  ; ++ir)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
              irzp =  particle.n_zgrid * particle.n_pgrid * ir +  particle.n_pgrid * iz + ip;
	   
	      particle.cr_density.d2[ir][iz].s[ip] =      particle_cr_density_arr[irzp];                    
	      
	    
	    }//ip
	   }//iz
	  }//ir



    //----------------------------------------------------------

			//Gulli20070810
    
    alpha1_r.delete_array();
    alpha1_z.delete_array();
    alpha1_p.delete_array();
    alpha2_r.delete_array();
    alpha2_z.delete_array();
    alpha2_p.delete_array();
    alpha3_r.delete_array();
    alpha3_z.delete_array();
    alpha3_p.delete_array();
    
    Nr1_.delete_array();
    Nz1_.delete_array();
    Np1_.delete_array();
    Nr2_.delete_array();
    Nz2_.delete_array();
    Np2_.delete_array();
    Nr3_.delete_array();
    Nz3_.delete_array();
    Np3_.delete_array();
    
    total_source_function.delete_array();
    Work.delete_array();


   if(galdef.solution_method==21)
   {
    delete[]      particle_cr_density_arr;                                                      //AWS20110204
    delete[]    total_source_function_arr;                                                      //AWS20110204
    delete[]                 alpha1_r_arr;                                                      //AWS20110204
    delete[]                 alpha1_z_arr;                                                      //AWS20110204
    delete[]                 alpha1_p_arr;                                                      //AWS20110204
    delete[]                 alpha2_r_arr;                                                      //AWS20110204
    delete[]                 alpha2_z_arr;                                                      //AWS20110204
    delete[]                 alpha2_p_arr;                                                      //AWS20110204
    delete[]                 alpha3_r_arr;                                                      //AWS20110204
    delete[]                 alpha3_z_arr;                                                      //AWS20110204
    delete[]                 alpha3_p_arr;                                                      //AWS20110204
    delete[]                     Work_arr;                                                      //AWS20110204
   }
 
  }   //2D case
  











  ////////////////////////////////////     3D     //////////////////////////////////////













  if (3 == particle.n_spatial_dimensions) {             // ==== 3D ====
   
    int ix,iy,iz,ip;
    int ixyzp,iixyzp;          //AWS20110202

    //Gulli20070810
    propel_arrays_initialize = 0;
    
    ostringstream aBuf;
    aBuf << "Generating alpha for 3D";
    INFO(aBuf.str());


    alpha1_x.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha1_y.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha1_z.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha1_p.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha2_x.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha2_y.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha2_z.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha2_p.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha3_x.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha3_y.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha3_z.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha3_p.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    
    Nx1_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Ny1_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Nz1_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Np1_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Nx2_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Ny2_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Nz2_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Np2_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Nx3_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Ny3_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Nz3_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Np3_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    
    total_source_function.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Work.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    
    if(galdef.solution_method==21)
    {
    particle_cr_density_arr  =new double[particle.n_xgrid * particle.n_ygrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110202
    total_source_function_arr=new double[particle.n_xgrid * particle.n_ygrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110202
    alpha1_x_arr             =new double[particle.n_xgrid * particle.n_ygrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110202
    alpha1_y_arr             =new double[particle.n_xgrid * particle.n_ygrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110202
    alpha1_z_arr             =new double[particle.n_xgrid * particle.n_ygrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110202
    alpha1_p_arr             =new double[particle.n_xgrid * particle.n_ygrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110202
    alpha2_x_arr             =new double[particle.n_xgrid * particle.n_ygrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110202
    alpha2_y_arr             =new double[particle.n_xgrid * particle.n_ygrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110202
    alpha2_z_arr             =new double[particle.n_xgrid * particle.n_ygrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110202
    alpha2_p_arr             =new double[particle.n_xgrid * particle.n_ygrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110202
    alpha3_x_arr             =new double[particle.n_xgrid * particle.n_ygrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110202
    alpha3_y_arr             =new double[particle.n_xgrid * particle.n_ygrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110202
    alpha3_z_arr             =new double[particle.n_xgrid * particle.n_ygrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110202
    alpha3_p_arr             =new double[particle.n_xgrid * particle.n_ygrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110202    
    Work_arr                 =new double[particle.n_xgrid * particle.n_ygrid * particle.n_zgrid * particle.n_pgrid]; //AWS20110202
    }


      // linear arrays for more efficiency in method 2                                  AWS20110202
    if(galdef.solution_method==21)
    {
      ixyzp=0;                                                                        //AWS20110202
      for (ix = 0; ix < particle.n_xgrid  ; ++ix)                                     //AWS20110202
      for (iy = 0; iy < particle.n_ygrid  ; ++iy)                                     //AWS20110202
      for (iz = 0; iz < particle.n_zgrid  ; ++iz)                                     //AWS20110202          
      for (ip = 0; ip < particle.n_pgrid  ; ++ip)                                     //AWS20110202                
      {                                                                               //AWS20110202

         particle_cr_density_arr[ixyzp] =   particle.cr_density.d3[ix][iy][iz].s[ip]; //AWS20110202

         ixyzp++;                                                                     //AWS20110202
      }                                                                               //AWS20110202
    }// if solution_method==21

    
    alpha1_x = alpha1_y = alpha1_z = alpha1_p = 0;
    alpha2_x = alpha2_y = alpha2_z = alpha2_p = 0;
    alpha3_x = alpha3_y = alpha3_z = alpha3_p = 0;
    
    // DIFFUSION term / Dxx = const in space 

    

    for (ix = 0; ix < particle.n_xgrid; ++ix) {
      
      for (iy = 0; iy < particle.n_ygrid; ++iy) {
	
	for (iz = 0; iz < particle.n_zgrid; ++iz) {

	  for (ip = 0; ip < particle.n_pgrid; ++ip) {

	    alpha1_x.d3[ix][iy][iz].s[ip]
	      = particle.Dxx.d3[ix][iy][iz].s[ip]*pow(particle.dx, -2.);
	    
	    alpha1_y.d3[ix][iy][iz].s[ip]
	      = particle.Dxx.d3[ix][iy][iz].s[ip]*pow(particle.dy, -2.);
              
     //	    alpha1_z.d3[ix][iy][iz].s[ip]                               //AWS20131017
     //	      = particle.Dxx.d3[ix][iy][iz].s[ip]*pow(particle.dz, -2.);//AWS20131017

	    alpha1_z.d3[ix][iy][iz].s[ip]                               //AWS20131017
	      = particle.Dzz.d3[ix][iy][iz].s[ip]*pow(particle.dz, -2.);//AWS20131017
              
	  }   //ip
          
	}   //iz
        
      }   //iy
      
    }   //ix

    alpha2_x = alpha1_x;  
    alpha2_x *= 2.;
    
    alpha3_x = alpha1_x;
    
    alpha2_y = alpha1_y;  
    alpha2_y *= 2.;
    
    alpha3_y = alpha1_y;
    
    alpha2_z = alpha1_z;  
    alpha2_z *= 2.;
    
    alpha3_z = alpha1_z;

    alpha1_x *= factor;   
    alpha1_y *= factor;   
    alpha1_z *= factor;
    
    alpha2_x *= factor;   
    alpha2_y *= factor;   
    alpha2_z *= factor; 
    
    alpha3_x *= factor;   
    alpha3_y *= factor;   
    alpha3_z *= factor;

// CONVECTION                                                      AWS20010330

    if (galdef.convection==1) {
     
      for (iz = 0; iz < particle.n_zgrid; ++iz) { // numerator = abs convection velocity in cm s^-1
		


	/* replaced by v_conv from array AWS20131008
	double az1 = (particle.z[iz] > 0.) ? 
                     (galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz-1]))*1.e5/particle.dz/kpc2cm: 0.; //AWS20110330 removed dt

	double az2 = (galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz]  ))*1.e5/particle.dz/kpc2cm;     //AWS20110330 removed dt
	
	double az3 = (particle.z[iz] < 0.) ? 
                     (galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz+1]))*1.e5/particle.dz/kpc2cm: 0.;//AWS20110330 removed dt
	*/


	// should be exactly the same for v0_conv=0

	double az1 = (particle.z[iz] > 0.) ? fabs(particle.v_conv.d3[0][0][iz-1].s[0])*1.e5/particle.dz/kpc2cm: 0.; //AWS20131008

	double az2 =                         fabs(particle.v_conv.d3[0][0][iz  ].s[0])*1.e5/particle.dz/kpc2cm;     //AWS20131008
	
	double az3 = (particle.z[iz] < 0.) ? fabs(particle.v_conv.d3[0][0][iz+1].s[0])*1.e5/particle.dz/kpc2cm: 0.; //AWS20131008


	if (fabs(particle.z[iz]) < 0.5*particle.dz) 
	  az1 = az2 = az3 = 0.;
	
        az1 *= dt; //AWS20110330
        az2 *= dt; //AWS20110330
        az3 *= dt; //AWS20110330



	for (ip = 0; ip < particle.n_pgrid-1; ++ip) {

	 /* replaced with dvdz array AWS20131008
	  double ap2 = particle.p[ip]/3.*galdef.dvdz_conv/(particle.p[ip+1] - particle.p[ip])*1.e5/kpc2cm;
	  
	  double ap3 = particle.p[ip+1]/3.*galdef.dvdz_conv/(particle.p[ip+1] - particle.p[ip])*1.e5/kpc2cm;
	 */

	  // keep same scheme (use ix=iy=0) for consistency with original, but later can generalize to convection varying with x,y


          ix=0;iy=0;                                                                                                                      //AWS20131008
	  const double ap2 = particle.p[ip]  /3.*particle.dvdz_conv.d3[ix][iy][iz].s[ip] / (particle.p[ip+1]-particle.p[ip])*1.e5/kpc2cm; //AWS20131008
	  
	  const double ap3 = particle.p[ip+1]/3.*particle.dvdz_conv.d3[ix][iy][iz].s[ip]/  (particle.p[ip+1]-particle.p[ip])*1.e5/kpc2cm; //AWS20131008

	  for (ix = 0; ix < particle.n_xgrid; ++ix)
	    for (iy = 0; iy < particle.n_ygrid; ++iy) {

	      alpha1_z.d3[ix][iy][iz].s[ip] += az1;                                      
	      alpha2_z.d3[ix][iy][iz].s[ip] += az2;  
	      alpha2_p.d3[ix][iy][iz].s[ip] += ap2;
	
	      alpha3_z.d3[ix][iy][iz].s[ip] += az3;  
	      alpha3_p.d3[ix][iy][iz].s[ip] += ap3;
	      
	    }   //iy ix
            
	}   //ip
       
	// NB not set for last ip grid point (see 2D): to do ! AWS20131008 
      }   //iz
      
    }

    // Radial CONVECTION                                                      ECC20150715
    if (galdef.convection==4) {
     
     for (ix = 0; ix < particle.n_xgrid; ++ix) { // numerator = abs convection velocity in cm s^-1
      for (iy = 0; iy < particle.n_ygrid; ++iy) { // numerator = abs convection velocity in cm s^-1
        for (iz = 0; iz < particle.n_zgrid; ++iz) { // numerator = abs convection velocity in cm s^-1
    
          double ax1 = (particle.x[ix] > 0.) ? fabs(particle.v_conv_x.d3[ix-1][iy][iz].s[0])*1.e5/particle.dx/kpc2cm: 0.; 
          double ax2 =                         fabs(particle.v_conv_x.d3[ix  ][iy][iz].s[0])*1.e5/particle.dx/kpc2cm;     
          double ax3 = (particle.x[ix] < 0.) ? fabs(particle.v_conv_x.d3[ix+1][iy][iz].s[0])*1.e5/particle.dx/kpc2cm: 0.; 

          double ay1 = (particle.y[iy] > 0.) ? fabs(particle.v_conv_y.d3[ix][iy-1][iz].s[0])*1.e5/particle.dy/kpc2cm: 0.; 
          double ay2 =                         fabs(particle.v_conv_y.d3[ix][iy  ][iz].s[0])*1.e5/particle.dy/kpc2cm;     
          double ay3 = (particle.y[iy] < 0.) ? fabs(particle.v_conv_y.d3[ix][iy+1][iz].s[0])*1.e5/particle.dy/kpc2cm: 0.; 

          double az1 = (particle.z[iz] > 0.) ? fabs(particle.v_conv_z.d3[ix][iy][iz-1].s[0])*1.e5/particle.dz/kpc2cm: 0.; 
          double az2 =                         fabs(particle.v_conv_z.d3[ix][iy][iz  ].s[0])*1.e5/particle.dz/kpc2cm;     
          double az3 = (particle.z[iz] < 0.) ? fabs(particle.v_conv_z.d3[ix][iy][iz+1].s[0])*1.e5/particle.dz/kpc2cm: 0.; 

          // Galactic center is a special case. 
          if ( particle.z[iz]==0 && particle.y[iy]==0 && particle.x[ix]==0) {
            ax1 = fabs(particle.v_conv_x.d3[ix-1][iy][iz].s[0])*1.e5/particle.dx/kpc2cm;  
            ax2 = fabs(particle.v_conv_x.d3[ix][iy][iz].s[0])*1.e5/particle.dx/kpc2cm; 
            ax3 = fabs(particle.v_conv_x.d3[ix+1][iy][iz].s[0])*1.e5/particle.dx/kpc2cm; 

            ay1 = fabs(particle.v_conv_x.d3[ix][iy-1][iz].s[0])*1.e5/particle.dy/kpc2cm;  
            ay2 = fabs(particle.v_conv_x.d3[ix][iy][iz].s[0])*1.e5/particle.dy/kpc2cm; 
            ay3 = fabs(particle.v_conv_x.d3[ix][iy+1][iz].s[0])*1.e5/particle.dy/kpc2cm; 

            az1 = fabs(particle.v_conv_x.d3[ix][iy][iz-1].s[0])*1.e5/particle.dz/kpc2cm;  
            az2 = fabs(particle.v_conv_x.d3[ix][iy][iz].s[0])*1.e5/particle.dz/kpc2cm; 
            az3 = fabs(particle.v_conv_x.d3[ix][iy][iz+1].s[0])*1.e5/particle.dz/kpc2cm; 
          }

          ax1 *= dt;
          ax2 *= dt; 
          ax3 *= dt; 

          ay1 *= dt; 
          ay2 *= dt; 
          ay3 *= dt; 

          az1 *= dt; 
          az2 *= dt; 
          az3 *= dt; 
         
          if (isnan(ax2)){
          cout  << " alpha1_x, alpha2_x, alpha3_x before ----- ax1,ax2,ax3:" << alpha1_x.d3[ix][iy][iz].s[ip] << " " << alpha2_x.d3[ix][iy][iz].s[ip]<< " " << alpha3_x.d3[ix][iy][iz].s[ip] << " " <<ax1 << " " <<ax2 << " " <<ax3 << endl;
          cout << "x, y, z: " << particle.x[ix] << " " << particle.y[iy] << " " << particle.z[iz] <<  endl; 
          cout << "particle.v_convx: " <<  fabs(particle.v_conv_x.d3[ix][iy][iz].s[0]) << endl;


          }


          for (ip = 0; ip < particle.n_pgrid-1; ++ip) {
            alpha1_x.d3[ix][iy][iz].s[ip] += ax1;                                      
            alpha1_y.d3[ix][iy][iz].s[ip] += ay1;                                      
            alpha1_z.d3[ix][iy][iz].s[ip] += az1;     

            alpha2_x.d3[ix][iy][iz].s[ip] += ax2;                                      
            alpha2_y.d3[ix][iy][iz].s[ip] += ay2;                                      
            alpha2_z.d3[ix][iy][iz].s[ip] += az2;

            alpha3_x.d3[ix][iy][iz].s[ip] += ax3;                                      
            alpha3_y.d3[ix][iy][iz].s[ip] += ay3;                                      
            alpha3_z.d3[ix][iy][iz].s[ip] += az3;                                                                            



          // //======================================================================
          // //debug ---------------------

          // double vv = 100; 
          // if (particle.z[iz] < 0){
          //   vv = -100;
          // }
            
          // double az1_b = (particle.z[iz] > 0.) ? vv*1.e5/particle.dz/kpc2cm: 0.; //AWS20131008
          // double az2_b =                         vv*1.e5/particle.dz/kpc2cm;     //AWS20131008
          // double az3_b = (particle.z[iz] < 0.) ? vv*1.e5/particle.dz/kpc2cm: 0.; //AWS20131008
          
          // az1_b *= dt; 
          // az2_b *= dt; 
          // az3_b *= dt; 

          // cout  << " az1_b, ... , az1 :" << az1_b << " " << az2_b<< " " << az3_b << " " <<az1 << " " <<az2 << " " <<az3 << endl;
          // //debug ---------------------

          } // ip 
      }   //iz
    } // iy
  } // ix
      
    }





// DIFFUSIVE REACCELERATION   IMOS20020329

    if (galdef.diff_reacc>0) { //IMOS20060330
     
      
      for (ix = 0; ix < particle.n_xgrid; ++ix) {
	
	for (iy = 0; iy < particle.n_ygrid; ++iy) {
	  
	  for (iz = 0; iz < particle.n_zgrid; ++iz) {

	    for (ip = 1; ip < particle.n_pgrid-1; ++ip) {
	      // alternative scheme #2, most detailed
	      alpha1_p.d3[ix][iy][iz].s[ip] += 
		(-(particle.Dpp.d3[ix][iy][iz].s[ip] - 
		   particle.Dpp.d3[ix][iy][iz].s[ip-1])
		 /(particle.p[ip] - particle.p[ip-1])
		 + 2.*particle.Dpp.d3[ix][iy][iz].s[ip]
		 /(particle.p[ip+1]-particle.p[ip-1])
		 + 2.*particle.Dpp.d3[ix][iy][iz].s[ip-1]
		 /particle.p[ip-1]) 
		/(particle.p[ip] - particle.p[ip-1]);

	      alpha2_p.d3[ix][iy][iz].s[ip] +=
		-(particle.Dpp.d3[ix][iy][iz].s[ip] - 
		  particle.Dpp.d3[ix][iy][iz].s[ip-1])
		/pow(particle.p[ip] - particle.p[ip-1], 2.)
		+ 2.*particle.Dpp.d3[ix][iy][iz].s[ip]
		/(particle.p[ip+1] - particle.p[ip-1])
		*(1./(particle.p[ip+1] - particle.p[ip])
		  + 1./(particle.p[ip] - particle.p[ip-1]))
		+ 2.*particle.Dpp.d3[ix][iy][iz].s[ip]
		/(particle.p[ip] - particle.p[ip-1])
		/particle.p[ip]; 
                  
	      alpha3_p.d3[ix][iy][iz].s[ip] +=
		2.*particle.Dpp.d3[ix][iy][iz].s[ip]
		/(particle.p[ip+1] - particle.p[ip-1])
		/(particle.p[ip+1] - particle.p[ip]);

	    
	    }   //ip
	  
	  }   //iz
          
	}   //iy
        
      }   //ix
      
    }   // diffusive reacceleration

    // MOMENTUM LOSSES
    
    if(galdef.momentum_losses) {
     
      for (ix = 0; ix < particle.n_xgrid; ++ix) {

	for (iy = 0; iy < particle.n_ygrid; ++iy) {

	  for (iz = 0; iz < particle.n_zgrid; ++iz) {

	    for (ip = 1; ip < particle.n_pgrid-1; ++ip) {// bug fixed IMOS20020417
	      
	      alpha2_p.d3[ix][iy][iz].s[ip] += 
		particle.dpdt.d3[ix][iy][iz].s[ip]
		/(particle.p[ip+1] - particle.p[ip]);
                     
	      alpha3_p.d3[ix][iy][iz].s[ip] += 
		particle.dpdt.d3[ix][iy][iz].s[ip+1]
		/(particle.p[ip+1]-particle.p[ip]);
                
	    }   //ip
            
	  }   //iz
          
	}   //iy
        
      }   //ix 
      
    }   //momentum losses

    // fill in values at extremes of p grid

    
    for (ix = 0; ix < particle.n_xgrid; ++ix) {
      
      for (iy = 0; iy < particle.n_ygrid; ++iy) {

	for (iz = 0; iz < particle.n_zgrid; ++iz) {

	  alpha1_p.d3[ix][iy][iz].s[0] = alpha1_p.d3[ix][iy][iz].s[1];
	  alpha2_p.d3[ix][iy][iz].s[0] = alpha2_p.d3[ix][iy][iz].s[1];
	  alpha3_p.d3[ix][iy][iz].s[0] = alpha3_p.d3[ix][iy][iz].s[1];
	  
	  alpha1_p.d3[ix][iy][iz].s[particle.n_pgrid-1] = alpha1_p.d3[ix][iy][iz].s[particle.n_pgrid-2];
	  alpha2_p.d3[ix][iy][iz].s[particle.n_pgrid-1] = alpha2_p.d3[ix][iy][iz].s[particle.n_pgrid-2];
	  alpha3_p.d3[ix][iy][iz].s[particle.n_pgrid-1] = alpha3_p.d3[ix][iy][iz].s[particle.n_pgrid-2];
          
	}   //iz
      
      }   //iy
      
    }   //ix 

    alpha1_p *= dt;
    alpha2_p *= dt;
    alpha3_p *= dt;

    double f_use = galdef.prop_x + galdef.prop_y + galdef.prop_z + galdef.prop_p;
    


    /* changing to method-independent alpha definition without f_use AWS20110120
    // FRAGMENTATION

    if (galdef.fragmentation) {
       
      Work = particle.fragment;        // use Work to avoid memory leak
      Work *= (dt/f_use);
      alpha2_x += Work;                 //  fragment used f_use times
      alpha2_y += Work;  
      alpha2_z += Work;                         
      alpha2_p += Work;

    
    }

    // RADIOACTIVE DECAY

    if (particle.t_half != 0.0 && galdef.radioactive_decay) {
         
      ostringstream rBuf;
      rBuf << "Radioactive decay of "<<particle.name<<" with half_life "
	   <<particle.t_half;
      INFO(rBuf.str());

      Work = particle.decay;           // use Work to avoid memory leak
      Work *= (dt/f_use);
      alpha2_x += Work;                 //  fragment used f_use times
      alpha2_y += Work; 
      alpha2_z += Work;                         
      alpha2_p += Work; 


    }
    */


    //AWS20110120
    //method-independent alpha definitions without f_use, using alpha_2 only once AWS20110220
    // assign to exactly one of the dimensions, preferring momentum if present, otherwise x or y or z

    // FRAGMENTATION
    
    if (galdef.fragmentation) {
       
      Work  = particle.fragment;         // use Work to avoid memory leak
      Work *= dt;

      if      (galdef.prop_p)alpha2_p += Work; 
      else if (galdef.prop_x)alpha2_x += Work;  
      else if (galdef.prop_y)alpha2_y += Work;                  
      else if (galdef.prop_z)alpha2_z += Work;                         
                          
 
    }

    // DECAY
    
    if (particle.t_half != 0 && galdef.radioactive_decay) {

      ostringstream rBuf;
      rBuf << "Radioactive decay of " << particle.name << " with half life " << particle.t_half;
      INFO(rBuf.str());
 
      Work  = particle.decay;            // use Work to avoid memory leak
      Work *= dt;

      if      (galdef.prop_p)alpha2_p += Work; 
      else if (galdef.prop_x)alpha2_x += Work;  
      else if (galdef.prop_y)alpha2_y += Work;                  
      else if (galdef.prop_z)alpha2_z += Work; 
 
    }

    buf.str("");
    buf<<"alpha1_x:"<< alpha1_x.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha2_x:"<< alpha2_x.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha3_x:"<< alpha3_x.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_y:"<< alpha1_y.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha2_y:"<< alpha2_y.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha3_y:"<< alpha3_y.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_z:"<< alpha1_z.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha2_z:"<< alpha2_z.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha3_z:"<< alpha3_z.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_p:"<< alpha1_p.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha2_p:"<< alpha2_p.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha3_p:"<< alpha3_p.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());

    if (galdef.verbose == 10) {
       
      cout<<"alpha1_x:"<<endl;   alpha1_x.print();
      cout<<"alpha2_x:"<<endl;   alpha2_x.print();
      cout<<"alpha3_x:"<<endl;   alpha3_x.print();
      cout<<"alpha1_y:"<<endl;   alpha1_y.print();
      cout<<"alpha2_y:"<<endl;   alpha2_y.print();
      cout<<"alpha3_y:"<<endl;   alpha3_y.print();
      cout<<"alpha1_z:"<<endl;   alpha1_z.print();
      cout<<"alpha2_z:"<<endl;   alpha2_z.print();
      cout<<"alpha3_z:"<<endl;   alpha3_z.print();
      cout<<"alpha1_p:"<<endl;   alpha1_p.print();
      cout<<"alpha2_p:"<<endl;   alpha2_p.print();
      cout<<"alpha3_p:"<<endl;   alpha3_p.print();

    }
    
    ostringstream tBuf;

    tBuf <<"galdef_ID "<<galdef.galdef_ID<<" start timestep dt="<<dt<<"  "
         <<particle.name;
    DEBUGLOG(tBuf.str());

    int timestep_mode;
    int timestep_mode2done = 0;
    
    for (timestep_mode = 1; timestep_mode <= 2; ++timestep_mode) {
      
       buf.str("");
      buf<<"    timestep_mode="<<timestep_mode;
      INFO(buf.str());
      t = 0;
         
      while ((1 == start) || 
	     (dt > end_timestep_sec/galdef.timestep_factor*0.9999) || 
	     (2 == timestep_mode && 0 == timestep_mode2done)) {
         
	if (0 == start && 1 == timestep_mode) {

	  dt *= galdef.timestep_factor;
          
	  ostringstream nTBuf;
   
	  nTBuf<<" galdef_ID "<<galdef.galdef_ID<<" new timestep dt="<<dt<<"  "
	      <<particle.name;
	  DEBUGLOG(nTBuf.str());
	  
	  alpha1_x *= galdef.timestep_factor;
	  alpha2_x *= galdef.timestep_factor;
	  alpha3_x *= galdef.timestep_factor;
	  alpha1_y *= galdef.timestep_factor;
	  alpha2_y *= galdef.timestep_factor;
	  alpha3_y *= galdef.timestep_factor;
	  alpha1_z *= galdef.timestep_factor;
	  alpha2_z *= galdef.timestep_factor;
	  alpha3_z *= galdef.timestep_factor;
	  alpha1_p *= galdef.timestep_factor;
	  alpha2_p *= galdef.timestep_factor;
	  alpha3_p *= galdef.timestep_factor;
	
	}
	
	start = 0;
	
	

	
	Nx1_ = alpha1_x;   Nx1_*= -0.5;
	Nx2_ = alpha2_x;   Nx2_*=  0.5;   Nx2_+= 1.0;   
	Nx3_ = alpha3_x;   Nx3_*= -0.5;
	Ny1_ = alpha1_y;   Ny1_*= -0.5;
	Ny2_ = alpha2_y;   Ny2_*=  0.5;   Ny2_+= 1.0;   
	Ny3_ = alpha3_y;   Ny3_*= -0.5;
	Nz1_ = alpha1_z;   Nz1_*= -0.5;
	Nz2_ = alpha2_z;   Nz2_*=  0.5;   Nz2_+= 1.0;   
	Nz3_ = alpha3_z;   Nz3_*= -0.5;
	Np1_ = alpha1_p;   Np1_*= -0.5;
	Np2_ = alpha2_p;   Np2_*=  0.5;   Np2_+= 1.0;   
	Np3_ = alpha3_p;   Np3_*= -0.5;
	
	buf.str("");
	buf<<" Nx,Ny,Nz,Np recomputed ";
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Nx1_    :"<< Nx1_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Nx2_    :"<< Nx2_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Nx3_    :"<< Nx3_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Ny1_    :"<< Ny1_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
       	buf<<"Ny2_    :"<< Ny2_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Ny3_    :"<< Ny3_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Nz1_    :"<< Nz1_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Nz2_    :"<< Nz2_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Nz3_    :"<< Nz3_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Np1_    :"<< Np1_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Np2_    :"<< Np2_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Np3_    :"<< Np3_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	
      
	total_source_function  = particle  .primary_source_function;
	total_source_function += particle.secondary_source_function;
        total_source_function *= dt;                                            //AWS20110131 for consistency with alpha



      if(galdef.solution_method==21)
      {
      // valarrays for more efficiency in method 2                                      AWS20110202
      ixyzp=0;                                                                        //AWS20110202
      for (ix = 0; ix < particle.n_xgrid  ; ++ix)                                     //AWS20110202
      for (iy = 0; iy < particle.n_ygrid  ; ++iy)                                     //AWS20110202   
      for (iz = 0; iz < particle.n_zgrid  ; ++iz)                                     //AWS20110202          
      for (ip = 0; ip < particle.n_pgrid  ; ++ip)                                     //AWS20110202                
      {                                                                               //AWS20110202

       total_source_function_arr[ixyzp] = total_source_function.d3[ix][iy][iz].s[ip]; //AWS20110202

                    alpha1_x_arr[ixyzp] =              alpha1_x.d3[ix][iy][iz].s[ip]; //AWS20110202
                    alpha1_y_arr[ixyzp] =              alpha1_y.d3[ix][iy][iz].s[ip]; //AWS20110202
                    alpha1_z_arr[ixyzp] =              alpha1_z.d3[ix][iy][iz].s[ip]; //AWS20110202
                    alpha1_p_arr[ixyzp] =              alpha1_p.d3[ix][iy][iz].s[ip]; //AWS20110202
                    alpha2_x_arr[ixyzp] =              alpha2_x.d3[ix][iy][iz].s[ip]; //AWS20110202
                    alpha2_y_arr[ixyzp] =              alpha2_y.d3[ix][iy][iz].s[ip]; //AWS20110202
                    alpha2_z_arr[ixyzp] =              alpha2_z.d3[ix][iy][iz].s[ip]; //AWS20110202
                    alpha2_p_arr[ixyzp] =              alpha2_p.d3[ix][iy][iz].s[ip]; //AWS20110202
                    alpha3_x_arr[ixyzp] =              alpha3_x.d3[ix][iy][iz].s[ip]; //AWS20110202
                    alpha3_y_arr[ixyzp] =              alpha3_y.d3[ix][iy][iz].s[ip]; //AWS20110202
                    alpha3_z_arr[ixyzp] =              alpha3_z.d3[ix][iy][iz].s[ip]; //AWS20110202
                    alpha3_p_arr[ixyzp] =              alpha3_p.d3[ix][iy][iz].s[ip]; //AWS20110202

       ixyzp++;                                                                       //AWS20110202
      }                                                                               //AWS20110202
      }// if solution_method==21
       

     
	if (0 == total_source_function.max()) 
	  return 0;
	
	int nrept;
	if (1 == timestep_mode) 
	  nrept = galdef.timestep_repeat;

	if (2 == timestep_mode) 
	  nrept = galdef.timestep_repeat2;
	
	


	if (1 == timestep_mode)
        {
          method = 1;  // fixed to ensure stability
	 
	  buf.str(""); buf<<"timestep_mode 1: using solution method "<<method<< " timestep="<<dt/year2sec<< " yr";
          INFO(buf);       
        }

	if (2 == timestep_mode)
        {
	 
          if(galdef.solution_method== 1){method=1;         }    //AWS20110201
	  if(galdef.solution_method== 2){method=2; turbo=0;}    //AWS20110201
          if(galdef.solution_method==21){method=2; turbo=1;}    //AWS20110201

          if(method==1){buf.str(""); buf<<"timestep_mode 2: using  method "<<method<<" Crank-Nicolson"
                                                          << " timestep="<<dt/year2sec<< " yr";}
	  if(method==2){buf.str(""); buf<<"timestep_mode 2: using  method "<<method<<" fully time explicit, "
		                        <<" turbo="<<turbo<< " timestep="<<dt/year2sec<< " yr" ;}

          INFO(buf);   


	  // following required to transfer from method=1 in timestep_mode 1
         if(galdef.solution_method==21)
         {
          ixyzp=0;                                                                        //AWS20110202

          for (ix = 0; ix < particle.n_xgrid  ; ++ix)                                     //AWS20110202
          for (iy = 0; iy < particle.n_ygrid  ; ++iy)                                     //AWS20110202
          for (iz = 0; iz < particle.n_zgrid  ; ++iz)                                     //AWS20110202          
          for (ip = 0; ip < particle.n_pgrid  ; ++ip)                                     //AWS20110202                
          {                                                                               //AWS20110202
           particle_cr_density_arr[ixyzp] =   particle.cr_density.d3[ix][iy][iz].s[ip];   //AWS20110202
           ixyzp++;                                                                       //AWS20110202
          }                                                                               //AWS20110202
         }// if solution_method==21

    
        }


        int  irept;                                 //AWS20111018
	for (irept = 1; irept <= nrept; ++irept) {  //AWS20111018

    
    // Modification by ECC.  If source_parameters[9]==0, then we shut off the 
    // primary source for all timesteps other than the first one.
    // i.e. a delta function injection. 	  
    

    if (irept>1 && galdef.source_parameters[9]==0 && 2==timestep_mode){
        if (irept%10==0) {cout << "Timestep:"<< irept <<" Injection Mode: " <<galdef.source_parameters[9]<<endl;}
        total_source_function=0;
        particle.primary_source_function=0;
        total_source_function += particle.secondary_source_function;
    }


	if (method == 1) // Crank-Nicolson method     AWS20111018
	{

	  if (1 == galdef.SNR_events && 2 == timestep_mode) {

	    t += dt;   
	    buf.str();
	    buf<<"timestep_mode 2: t="<<t;
	    INFO(buf.str());
	    source_SNR_event(particle, t/year2sec); // AWS20000807
	    total_source_function  = particle.primary_source_function;
	    total_source_function += particle.secondary_source_function;
            total_source_function *= dt;                                        //AWS20110131 for consistency with alpha

	    if(galdef.verbose>=1) {

	       buf.str("");
	      buf<<"propel:source_SNR_event: primary source function for particle " <<particle.name;
	      DEBUGLOG(buf.str());
	      particle.primary_source_function .print();
	    
	    }
            
	  }   // SNR_events

	  int n_xgrid_sym = particle.n_xgrid;
	  int n_ygrid_sym = particle.n_ygrid;
	  int n_zgrid_sym = particle.n_zgrid;

	  if (galdef.use_symmetry == 2) {
	  
	    n_xgrid_sym = particle.n_xgrid/2;
	    n_ygrid_sym = particle.n_ygrid/2;
	    n_zgrid_sym = particle.n_zgrid/2;
            
	  }
	  
	 
	  
	 

	    // X propagation
	    
 //Gulli20070821
#pragma omp parallel default(shared) private(iz,ix,iy,ip)
	    if (1 == galdef.prop_x) {



	      valarray<double> Nx0(0., particle.n_xgrid);
	      valarray<double> Nx1(0., particle.n_xgrid);
	      valarray<double> Nx2(0., particle.n_xgrid);
	      valarray<double> Nx3(0., particle.n_xgrid);
	      valarray<double> Rx (0., particle.n_xgrid);
	      
 //Gulli20070821
#pragma omp for schedule(runtime)
	      for (ip = 0; ip < particle.n_pgrid; ++ip) {
		
		for (iz = 0; iz < particle.n_zgrid; ++iz) {
		  
		  for (iy = 0; iy < particle.n_ygrid; ++iy) {

		    if (galdef.use_symmetry <= 1 || 
			(iy <= n_ygrid_sym && 
			 iz <= n_zgrid_sym)) {
                         
		      for (ix = 0; ix < particle.n_xgrid; ++ix) {

			Nx1[ix] = Nx1_.d3[ix][iy][iz].s[ip];
			Nx2[ix] = Nx2_.d3[ix][iy][iz].s[ip];
			Nx3[ix] = Nx3_.d3[ix][iy][iz].s[ip];
			Nx0[ix] = total_source_function.d3[ix][iy][iz].s[ip]   /f_use;//AWS20110131

			Nx0[ix] += (1. - 0.5*alpha2_x.d3[ix][iy][iz].s[ip])*particle.cr_density.d3[ix][iy][iz].s[ip];
                     
			if (ix > 0) 
			  Nx0[ix] += 0.5*alpha1_x.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix-1][iy][iz].s[ip];
                    
			if (1 == galdef.use_symmetry && 0 == ix) 
			  Nx0[ix] += 0.5*alpha1_x.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix+1][iy][iz].s[ip];
                                 
			if (ix < particle.n_xgrid-1)
                                 
			  Nx0[ix] += 0.5*alpha3_x.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix+1][iy][iz].s[ip];
                              
		      }
                      
		      if (1 != galdef.use_symmetry) 
			tridag(&Nx1[0], &Nx2[0], &Nx3[0], &Nx0[0], &Rx[0], particle.n_xgrid);  
		      
		      if (1 == galdef.use_symmetry)
			tridag_sym(&Nx1[0], &Nx2[0], &Nx3[0], &Nx0[0], &Rx[0], particle.n_xgrid); 

		      for (ix = 0; ix < particle.n_xgrid; ++ix)
			particle.cr_density.d3[ix][iy][iz].s[ip] = (Rx[ix]< 0. ? 0.: Rx[ix]); //IMOS20030217
            
		    }   //symmetry
                    
		  }   //iy
                  
		}   //iz
                
	      }   //ip

	      if (2 == galdef.use_symmetry) {
 //Gulli20070821
#pragma omp for schedule(runtime)
		for (ip = 0; ip < particle.n_pgrid; ++ip) {

		  for (iz = 0; iz <= n_zgrid_sym; ++iz) {

		    for (iy = 0; iy <= n_ygrid_sym; ++iy) {

		      for (ix = 0; ix < particle.n_xgrid; ++ix) {

			double value = particle.cr_density.d3[ix][iy][iz].s[ip];  //Gulli20070810
			particle.cr_density.d3[ix][particle.n_ygrid-1-iy][iz].s[ip] = value;
			
			particle.cr_density.d3[ix][iy][particle.n_zgrid-1-iz].s[ip] = value;
			
			particle.cr_density.d3[ix][particle.n_ygrid-1-iy][particle.n_zgrid-1-iz].s[ip] = value;
                         
		      }   //ix
		    
		    }   //iy
                    
		  }   //iz
                  
		}   //ip
                
	      }

	    }   // prop_x

	    // Y propagation

	    //Gulli20070821
#pragma omp parallel default(shared) private(iz,ix,iy,ip)
	    if (1 == galdef.prop_y) {



	      valarray<double> Ny0(0., particle.n_ygrid);
	      valarray<double> Ny1(0., particle.n_ygrid);
	      valarray<double> Ny2(0., particle.n_ygrid);
	      valarray<double> Ny3(0., particle.n_ygrid);
	      valarray<double> Ry(0., particle.n_ygrid);

	      //Gulli20070821
#pragma omp for schedule(runtime)
	      for (ip = 0; ip < particle.n_pgrid; ++ip) {
		
		for (iz = 0; iz < particle.n_zgrid; ++iz) {
		  
		  for (ix = 0; ix < particle.n_xgrid; ++ix) {

		    if (galdef.use_symmetry <= 1 || 
			(iz <= n_zgrid_sym && 
			 ix <= n_xgrid_sym)) {

		      for (iy = 0; iy < particle.n_ygrid; ++iy) {

			Ny1[iy] = Ny1_.d3[ix][iy][iz].s[ip];
			Ny2[iy] = Ny2_.d3[ix][iy][iz].s[ip];
			Ny3[iy] = Ny3_.d3[ix][iy][iz].s[ip];
			Ny0[iy] = total_source_function.d3[ix][iy][iz].s[ip]   /f_use;//AWS20110131
                        
			Ny0[iy] += (1. - 0.5*alpha2_y.d3[ix][iy][iz].s[ip])*particle.cr_density.d3[ix][iy][iz].s[ip];
                         
			if (iy > 0) 
			  Ny0[iy] += 0.5*alpha1_y.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix][iy-1][iz].s[ip];

			if (1 == galdef.use_symmetry && 0 == iy) 
			  Ny0[iy] += 0.5*alpha1_y.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix][iy+1][iz].s[ip];

			if (iy < particle.n_ygrid-1) 
			  Ny0[iy] += 0.5*alpha3_y.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix][iy+1][iz].s[ip];
                              
		      }

		      if (1 != galdef.use_symmetry) 
			tridag(&Ny1[0], &Ny2[0], &Ny3[0], &Ny0[0], &Ry[0], particle.n_ygrid);  
		      
		      if (1 == galdef.use_symmetry) 
			tridag_sym(&Ny1[0], &Ny2[0], &Ny3[0], &Ny0[0], &Ry[0], particle.n_ygrid);  
		      
		      for (iy = 0; iy < particle.n_ygrid; ++iy) 
			particle.cr_density.d3[ix][iy][iz].s[ip] = (Ry[iy] < 0. ? 0.: Ry[iy]); //IMOS20030217
             
		    }  //  symmetry

		  }  //  ix
                  
		}  //  iz
                
	      }  //  ip
	      
	      if (2 == galdef.use_symmetry) {
		//Gulli20070821
#pragma omp for schedule(runtime)
		for (ip = 0; ip < particle.n_pgrid; ++ip) {

		  for (iz = 0; iz <= n_zgrid_sym; ++iz) {
		    
		    for (ix = 0; ix <= n_xgrid_sym; ++ix) {

		      for (iy = 0; iy < particle.n_ygrid; ++iy) {

			double value =  particle.cr_density.d3[ix][iy][iz].s[ip];
			particle.cr_density.d3[particle.n_xgrid-1-ix][iy][iz].s[ip] = value;
			particle.cr_density.d3[ix][iy][particle.n_zgrid-1-iz].s[ip] = value;
			particle.cr_density.d3[particle.n_xgrid-1-ix][iy][particle.n_zgrid-1-iz].s[ip] = value;
                         
		      }  //  iy    
                      
		    }  //  ix     
                    
		  }  //  iz
                  
		}  //  ip
                
	      }  //  symmetry

	    }  //  prop_y

	    // Z propagation
	    
	    //Gulli20070821
#pragma omp parallel default(shared) private(iz,ix,iy,ip)
	    if (1 == galdef.prop_z) {
	      

	      
	      valarray<double> Nz0(0., particle.n_zgrid);
	      valarray<double> Nz1(0., particle.n_zgrid);
	      valarray<double> Nz2(0., particle.n_zgrid);
	      valarray<double> Nz3(0., particle.n_zgrid);
	      valarray<double> Rz(0., particle.n_zgrid);
	      
	      //Gulli20070821
#pragma omp for schedule(runtime)
	      for (ip = 0; ip < particle.n_pgrid; ++ip) {

		for (ix = 0; ix < particle.n_xgrid; ++ix) {

		  for (iy = 0; iy < particle.n_ygrid; ++iy) {

		    if (galdef.use_symmetry <= 1 || 
			(ix <= n_xgrid_sym && 
			 iy <= n_ygrid_sym)) {

		      for (iz = 0; iz < particle.n_zgrid; ++iz) {

			Nz1[iz] = Nz1_.d3[ix][iy][iz].s[ip];
			Nz2[iz] = Nz2_.d3[ix][iy][iz].s[ip];
			Nz3[iz] = Nz3_.d3[ix][iy][iz].s[ip];
			Nz0[iz] = total_source_function.d3[ix][iy][iz].s[ip]   /f_use;//AWS20110131

			Nz0[iz] += (1. - 0.5*alpha2_z.d3[ix][iy][iz].s[ip])*particle.cr_density.d3[ix][iy][iz].s[ip];
                            
			if (iz > 0)
			  Nz0[iz] += 0.5*alpha1_z.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix][iy][iz-1].s[ip];
			
			if (galdef.use_symmetry && iz==0)
			  Nz0[iz] += 0.5*alpha1_z.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix][iy][iz+1].s[ip];

			if (iz <particle.n_zgrid-1)
			  Nz0[iz] += 0.5*alpha3_z.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix][iy][iz+1].s[ip];

		      } 
// cout<<"ix iy ip "<<ix<<" "<<iy<<" "<<ip<<endl;
                      
		      if (1 != galdef.use_symmetry) 
			tridag(&Nz1[0], &Nz2[0], &Nz3[0], &Nz0[0], &Rz[0], particle.n_zgrid);  
		      
		      if (1 == galdef.use_symmetry) 
			tridag_sym(&Nz1[0], &Nz2[0], &Nz3[0], &Nz0[0], &Rz[0], particle.n_zgrid);  

		      for (iz = 0; iz < particle.n_zgrid; ++iz) 
			particle.cr_density.d3[ix][iy][iz].s[ip] = (Rz[iz] < 0. ? 0.: Rz[iz]); //IMOS20030217
             
		    }  //  symmetry

		  }  //  iy
                  
		}  //  ix
                
	      }  //  ip
	      
	      if (2 == galdef.use_symmetry) {
		//Gulli20070821
#pragma omp for schedule(runtime)
		for (ip = 0; ip < particle.n_pgrid; ++ip) {

		  for (ix = 0; ix <= n_xgrid_sym; ++ix) {

		    for (iy = 0; iy <= n_ygrid_sym; ++iy) {

		      for (iz = 0; iz < particle.n_zgrid; ++iz) {

			double value = particle.cr_density.d3[ix][iy][iz].s[ip];
			
			particle.cr_density.d3[particle.n_xgrid-1-ix][iy][iz].s[ip] = value;
			particle.cr_density.d3[ix][particle.n_ygrid-1-iy][iz].s[ip] = value;
			particle.cr_density.d3[particle.n_xgrid-1-ix][particle.n_ygrid-1-iy][iz].s[ip] = value;
                          
		      }  //  iz
                      
		    }  //  iy 
                    
		  }  //  ix
                  
		}  //  ip
                
	      }  //  symmetry

	    }  //  prop_z

	    // P propagation

	    //Gulli20070821
#pragma omp parallel default(shared) private(iz,ix,iy,ip)
	    if (1 == galdef.prop_p) { 


	      valarray<double> Np0(0., particle.n_pgrid);
	      valarray<double> Np1(0., particle.n_pgrid);
	      valarray<double> Np2(0., particle.n_pgrid);
	      valarray<double> Np3(0., particle.n_pgrid);
	      valarray<double> Rp(0., particle.n_pgrid);

	      //Gulli20070821
#pragma omp for schedule(runtime)
	      for (ix = 0; ix < particle.n_xgrid; ++ix) {

		for (iy = 0; iy < particle.n_ygrid; ++iy) {

		  for (iz = 0; iz < particle.n_zgrid; ++iz) {
		    
		    if (galdef.use_symmetry <= 1 || 
			(ix <= n_xgrid_sym && 
			 iy <= n_ygrid_sym && 
			 iz <= n_zgrid_sym)) {

		      for (ip = 0; ip < particle.n_pgrid; ++ip) {

			Np1[ip] = Np1_.d3[ix][iy][iz].s[ip];
			Np2[ip] = Np2_.d3[ix][iy][iz].s[ip];
			Np3[ip] = Np3_.d3[ix][iy][iz].s[ip];
			
			Np0[ip] = total_source_function.d3[ix][iy][iz].s[ip]   /f_use;//AWS20110131
                                 
			Np0[ip] += (1. - 0.5*alpha2_p.d3[ix][iy][iz].s[ip])*particle.cr_density.d3[ix][iy][iz].s[ip];
			
			if (ip > 0)
			  Np0[ip] += 0.5*alpha1_p.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix][iy][iz].s[ip-1];
			
			if (ip < particle.n_pgrid-1)
			  Np0[ip] += 0.5*alpha3_p.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix][iy][iz].s[ip+1];
                              
		      }
		      
		      tridag(&Np1[0], &Np2[0], &Np3[0], &Np0[0], &Rp[0], particle.n_pgrid);  
		      
		      for (ip = 0; ip < particle.n_pgrid; ++ip) 
			particle.cr_density.d3[ix][iy][iz].s[ip] = (Rp[ip] < 0. ? 0.: Rp[ip]); //IMOS20030217
		    
		    }  //  symmetry

		  }  //  iz
		
		}  //  iy
                
	      }  //  ix
	       
	      if (2 == galdef.use_symmetry) {

		//Gulli20070821
#pragma omp for schedule(runtime)
		for (ip = 0; ip < particle.n_pgrid; ++ip) {

		  for (ix = 0; ix <= n_xgrid_sym; ++ix) {

		    for (iy = 0; iy <= n_ygrid_sym; ++iy) {

		      for (iz = 0; iz <= n_zgrid_sym; ++iz) {

			double value = particle.cr_density.d3[ix][iy][iz].s[ip];
			
			particle.cr_density.d3[ix][particle.n_ygrid-1-iy][iz].s[ip] = value;
			particle.cr_density.d3[ix][iy][particle.n_zgrid-1-iz].s[ip] = value;
			particle.cr_density.d3[ix][particle.n_ygrid-1-iy][particle.n_zgrid-1-iz].s[ip] = value;
			particle.cr_density.d3[particle.n_xgrid-1-ix][iy][iz].s[ip] = value;
			particle.cr_density.d3[particle.n_xgrid-1-ix][iy][particle.n_zgrid-1-iz].s[ip] = value;
			particle.cr_density.d3[particle.n_xgrid-1-ix][particle.n_ygrid-1-iy][iz].s[ip] = value;
			particle.cr_density.d3[particle.n_xgrid-1-ix][particle.n_ygrid-1-iy][particle.n_zgrid-1-iz].s[ip] = value; 
		      
		      }  //  iz   
		    
		    }  //  iy      
                    
		  }  //  ix
                  
		}  //  ip
                
	      }  //  symmetry

               
	    }  //  prop_p
	    
	 

	  


	  // set BOUNDARY CONDITIONS: zero on all boundaries
	  if( galdef.spatial_bound_conds==1) //AWS20131028
	  {
 
	  //Gulli20070821
#pragma omp parallel for schedule(runtime) default(shared) private(iz,ix,iy,ip)
	  for (ip = 0; ip < particle.n_pgrid; ++ip) {
	    
	    for (ix = 0; ix < particle.n_xgrid; ++ix) {

	      for (iy = 0; iy < particle.n_ygrid; ++iy) {

		if (1 != galdef.use_symmetry) { 
	
		  iz = 0;  
		  particle.cr_density.d3[ix][iy][iz].s[ip] = 0; 
		
		}
		  
		iz = particle.n_zgrid-1;
		particle.cr_density.d3[ix][iy][iz].s[ip] = 0;

	      }
	     
	    }
	     
	    for (ix = 0; ix < particle.n_xgrid; ++ix) {
	      
	      for (iz = 0; iz < particle.n_zgrid; ++iz) {
		  
		if (1 != galdef.use_symmetry) { 
		  
		  iy = 0;  
		  particle.cr_density.d3[ix][iy][iz].s[ip] = 0; 

		}
		
		iy = particle.n_ygrid-1;
		particle.cr_density.d3[ix][iy][iz].s[ip] = 0;

	      }	      

	    }
	     
	    for (iz = 0; iz < particle.n_zgrid; ++iz) {
	      
	      for (iy = 0; iy < particle.n_ygrid; ++iy) {
		
		if (1 != galdef.use_symmetry) { 

		  ix = 0;  
		  particle.cr_density.d3[ix][iy][iz].s[ip] = 0; 

		}
		
		ix = particle.n_xgrid-1;
		particle.cr_density.d3[ix][iy][iz].s[ip] = 0;
		
	      }
	    
	    }

	  }  //  boundary conditions ip loop
	  
	 } // spatial_bound_conds==1

	 // end of updating step for Crank-Nicolson method
       } // method == 1
 






	if (method == 2) // explicit updating method AWS20110119
	{
	 

	  int nxyzp =  particle.n_xgrid *  particle.n_ygrid *  particle.n_zgrid *  particle.n_pgrid; //AWS20110202


	 if(turbo==0)                                                                               //AWS20110202
	 {

	  // group by boundaries for efficient optimization avoiding conditionals

#pragma omp parallel for schedule(runtime) private(ix,iy,iz,ip) default(shared)



	 for (ix = 0; ix < particle.n_xgrid  ; ++ix)
         {
	  for (iy = 0; iy < particle.n_ygrid  ; ++iy)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {

	      Work.d3[ix][iy][iz].s[ip]  =                               total_source_function.d3[ix][iy][iz].s[ip];//AWS20110131

	      Work.d3[ix][iy][iz].s[ip] -= alpha2_x.d3[ix][iy][iz].s[ip] * particle.cr_density.d3[ix][iy][iz].s[ip];
	      Work.d3[ix][iy][iz].s[ip] -= alpha2_y.d3[ix][iy][iz].s[ip] * particle.cr_density.d3[ix][iy][iz].s[ip];
	      Work.d3[ix][iy][iz].s[ip] -= alpha2_z.d3[ix][iy][iz].s[ip] * particle.cr_density.d3[ix][iy][iz].s[ip];  
	      Work.d3[ix][iy][iz].s[ip] -= alpha2_p.d3[ix][iy][iz].s[ip] * particle.cr_density.d3[ix][iy][iz].s[ip];
	      

	      /* removed to assist optimation/parallelization
              if(0==1)
	      {
              cout<<"method 2: dt="<<dt;
              cout<<"          total_source_function*dt="<<total_source_function.d3[ix][iy][iz].s[ip]*dt<<" " ;
              cout<<"          alpha1_r="<<alpha1_x.d3[ix][iy][iz].s[ip]<<" " ;
              cout<<"          alpha2_r="<<alpha2_x.d3[ix][iy][iz].s[ip]<<" " ;
              cout<<"          alpha3_r="<<alpha3_x.d3[ix][iy][iz].s[ip]<<" " ;
              cout<<"particle.cr_density="<<particle.cr_density.d3[ix][iy][iz].s[ip]<<" " ;
              cout<<"          Work=    "<<    Work.d3[ix][iy][iz].s[ip]<<endl;
	      }
	      */

	    }//ip
	   }//iz
	  }//iy
         }//ix



#pragma omp parallel for schedule(runtime) private(ix,iy,iz,ip) default(shared)
	 for (ix = 1; ix < particle.n_xgrid  ; ++ix)
         {
	  for (iy = 0; iy < particle.n_ygrid  ; ++iy)
          {
	    
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	      
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
	      Work.d3[ix][iy][iz].s[ip] += alpha1_x.d3[ix][iy][iz].s[ip] * particle.cr_density.d3[ix-1][iy][iz].s[ip];
	    }//ip
	   }//iz
	  }//iy
	 }//ix

#pragma omp parallel for schedule(runtime) private(ix,iy,iz,ip) default(shared)
	 for (ix = 0; ix < particle.n_xgrid-1; ++ix)
         {
	  for (iy = 0; iy < particle.n_ygrid  ; ++iy)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
	      Work.d3[ix][iy][iz].s[ip] += alpha3_x.d3[ix][iy][iz].s[ip] * particle.cr_density.d3[ix+1][iy][iz].s[ip];
	    }//ip
	   }//iz
	  }//iy
	 }//ix



#pragma omp parallel for schedule(runtime) private(ix,iy,iz,ip) default(shared)
	 for (ix = 0; ix < particle.n_xgrid  ; ++ix)
         {
	  for (iy = 1; iy < particle.n_ygrid  ; ++iy)
          {
	    
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	      
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
	      Work.d3[ix][iy][iz].s[ip] += alpha1_y.d3[ix][iy][iz].s[ip] * particle.cr_density.d3[ix][iy-1][iz].s[ip];
	    }//ip
	   }//iz
	  }//iy
	 }//ix

#pragma omp parallel for schedule(runtime) private(ix,iy,iz,ip) default(shared)
	 for (ix = 0; ix < particle.n_xgrid  ; ++ix)
         {
	  for (iy = 0; iy < particle.n_ygrid-1; ++iy)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
	      Work.d3[ix][iy][iz].s[ip] += alpha3_y.d3[ix][iy][iz].s[ip] * particle.cr_density.d3[ix][iy+1][iz].s[ip];
	    }//ip
	   }//iz
	  }//iy
	 }//ix


#pragma omp parallel for schedule(runtime) private(ix,iy,iz,ip) default(shared)
	 for (ix = 0; ix < particle.n_xgrid  ; ++ix)
         {
	  for (iy = 0; iy < particle.n_ygrid  ; ++iy)
          {
	    
	   for (iz = 1; iz < particle.n_zgrid  ; ++iz)
           {
	      
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
	      Work.d3[ix][iy][iz].s[ip] += alpha1_z.d3[ix][iy][iz].s[ip] * particle.cr_density.d3[ix][iy][iz-1].s[ip];
	    }//ip
	   }//iz
	  }//iy
	 }//ix

#pragma omp parallel for schedule(runtime) private(ix,iy,iz,ip) default(shared)
	 for (ix = 0; ix < particle.n_xgrid  ; ++ix)
         {
	  for (iy = 0; iy < particle.n_ygrid  ; ++iy)
          {
	   for (iz = 0; iz < particle.n_zgrid-1; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
	      Work.d3[ix][iy][iz].s[ip] += alpha3_z.d3[ix][iy][iz].s[ip] * particle.cr_density.d3[ix][iy][iz+1].s[ip];
	    }//ip
	   }//iz
	  }//iy
	 }//ix


#pragma omp parallel for schedule(runtime) private(ix,iy,iz,ip) default(shared)
	 for (ix = 0; ix < particle.n_xgrid  ; ++ix)
         {
	  for (iy = 0; iy < particle.n_ygrid  ; ++iy)
          {
	    
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	      
	    for (ip = 1; ip < particle.n_pgrid  ; ++ip)
            {
	      Work.d3[ix][iy][iz].s[ip] += alpha1_p.d3[ix][iy][iz].s[ip] * particle.cr_density.d3[ix][iy][iz].s[ip-1];
	    }//ip
	   }//iz
	  }//iy
	 }//ix

#pragma omp parallel for schedule(runtime) private(ix,iy,iz,ip) default(shared)
	 for (ix = 0; ix < particle.n_xgrid  ; ++ix)
         {
	  for (iy = 0; iy < particle.n_ygrid  ; ++iy)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid-1; ++ip)
            {
	      Work.d3[ix][iy][iz].s[ip] += alpha3_p.d3[ix][iy][iz].s[ip] * particle.cr_density.d3[ix][iy][iz].s[ip+1];
	    }//ip
	   }//iz
	  }//iy
	 }//ix
	 

	 //	  particle.cr_density += Work;  //AWS20110126 replaced by parallel version below

#pragma omp parallel for schedule(runtime) private(ix,iy,iz,ip) default(shared)
	 for (ix = 0; ix < particle.n_xgrid  ; ++ix)
         {
	  for (iy = 0; iy < particle.n_ygrid  ; ++iy)
          {
	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)
           {
	    for (ip = 0; ip < particle.n_pgrid  ; ++ip)
            {
	     particle.cr_density.d3[ix][iy][iz].s[ip] += Work.d3[ix][iy][iz].s[ip];
	    }//ip
	   }//iz
	  }//iy
	 }//ix



	  // apply spatial boundary conditions
	  if( galdef.spatial_bound_conds==1) //AWS20131028
	  { 

#pragma omp parallel for schedule(runtime) private(ix,iy,iz,ip) default(shared)
	  for (ip = 0; ip < particle.n_pgrid; ++ip) {
	    
	    for (ix = 0; ix < particle.n_xgrid; ++ix) {

	      for (iy = 0; iy < particle.n_ygrid; ++iy) {

		 
	
		  iz = 0;  
		  particle.cr_density.d3[ix][iy][iz].s[ip] = 0; 
		
		
		  
		iz = particle.n_zgrid-1;
		particle.cr_density.d3[ix][iy][iz].s[ip] = 0;

	      }
	     
	    }

     
	    for (ix = 0; ix < particle.n_xgrid; ++ix) {
	      
	      for (iz = 0; iz < particle.n_zgrid; ++iz) {
		  
	       
		  
		  iy = 0;  
		  particle.cr_density.d3[ix][iy][iz].s[ip] = 0; 

		
		
		iy = particle.n_ygrid-1;
		particle.cr_density.d3[ix][iy][iz].s[ip] = 0;

	      }	      

	    }
	     

	    for (iz = 0; iz < particle.n_zgrid; ++iz) {
	      
	      for (iy = 0; iy < particle.n_ygrid; ++iy) {
		
	        

		  ix = 0;  
		  particle.cr_density.d3[ix][iy][iz].s[ip] = 0; 

		
		
		ix = particle.n_xgrid-1;
		particle.cr_density.d3[ix][iy][iz].s[ip] = 0;
		
	      }
	    
	    }

	  }  // ip  

	  } // spatial_bound_conds==1




	} // turbo == 0                                                                                                                //AWS20110202





	  //----------------------------------------------------------------------------


	  if (turbo==1)                                                                                                                  //AWS20110202 
	 { 
	  
	   #pragma omp parallel for schedule(runtime)   default(shared)
	    for (ixyzp=0; ixyzp<nxyzp  ; ++ixyzp)
	    {
 

	      Work_arr[ixyzp]        =  total_source_function_arr[ixyzp]                                                                  //AWS20110202
                                                   -(alpha2_x_arr[ixyzp]+alpha2_y_arr[ixyzp]+alpha2_z_arr[ixyzp]+ alpha2_p_arr[ixyzp])    //AWS20110202
                                                  *                               particle_cr_density_arr[ixyzp];                         //AWS20110202
   
      
	    }//ixyzp                                                                                                                      //AWS20110202


	  // x direction -----------------------------------------------------------
             	  
#pragma omp parallel for schedule(runtime)  private(ix,iy,iz,ip,ixyzp,iixyzp) default(shared)       //AWS20110202
	   
	 for   (ix = 1; ix < particle.n_xgrid  ; ++ix)                                                //AWS20110202
         {
	  for  (iy = 0; iy < particle.n_ygrid  ; ++iy)                                               //AWS20110202
          {
  	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)                                              //AWS20110202
           {
             ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid * iz; //AWS20110202
	    iixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid * (ix-1)  + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid * iz; //AWS20110202
	    
            for (ip = 0; ip < particle.n_pgrid  ; ++ip)                                              //AWS20110202
            {
	     Work_arr[ixyzp]        += alpha1_x_arr[ixyzp]        * particle_cr_density_arr[iixyzp]; //AWS20110202

                      ixyzp++;                                                                       //AWS20110202
	             iixyzp++;                                                                       //AWS20110202

	    }//ip               	    
	   }//iz
	  }//iy
 	 }//ix

#pragma omp parallel for schedule(runtime)  private(ix,iy,iz,ip,ixyzp,iixyzp) default(shared)       //AWS20110202
	   
	 for   (ix = 0; ix < particle.n_xgrid-1; ++ix)                                                //AWS20110202
         {
	  for  (iy = 0; iy < particle.n_ygrid  ; ++iy)                                               //AWS20110202
          {
  	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)                                              //AWS20110202
           {
             ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid * iz; //AWS20110202
	    iixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid * (ix+1)  + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid * iz; //AWS20110202
	    
            for (ip = 0; ip < particle.n_pgrid  ; ++ip)                                              //AWS20110202
            {
	     Work_arr[ixyzp]        += alpha3_x_arr[ixyzp]        * particle_cr_density_arr[iixyzp]; //AWS20110202

                      ixyzp++;                                                                       //AWS20110202
	             iixyzp++;                                                                       //AWS20110202

	    }//ip               	    
	   }//iz
	  }//iy
 	 }//ix


	  // y direction -----------------------------------------------------------
             	  
#pragma omp parallel for schedule(runtime)  private(ix,iy,iz,ip,ixyzp,iixyzp) default(shared)       //AWS20110202
	   
	 for   (ix = 0; ix < particle.n_xgrid  ; ++ix)                                                //AWS20110202
         {
	  for  (iy = 1; iy < particle.n_ygrid  ; ++iy)                                               //AWS20110202
          {
  	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)                                              //AWS20110202
           {
             ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid * iz; //AWS20110202
	     iixyzp =  particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid * (iy-1) +  particle.n_pgrid * iz; //AWS20110202
	    
            for (ip = 0; ip < particle.n_pgrid  ; ++ip)                                              //AWS20110202
            {
	     Work_arr[ixyzp]        += alpha1_y_arr[ixyzp]        * particle_cr_density_arr[iixyzp]; //AWS20110202

                      ixyzp++;                                                                       //AWS20110202
	             iixyzp++;                                                                       //AWS20110202

	    }//ip               	    
	   }//iz
	  }//iy
 	 }//ix

#pragma omp parallel for schedule(runtime)  private(ix,iy,iz,ip,ixyzp,iixyzp) default(shared)       //AWS20110202
	   
	 for   (ix = 0; ix < particle.n_xgrid  ; ++ix)                                                //AWS20110202
         {
	  for  (iy = 0; iy < particle.n_ygrid-1; ++iy)                                               //AWS20110202
          {
  	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)                                              //AWS20110202
           {
             ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid * iz; //AWS20110202
	     iixyzp =  particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid * (iy+1) +  particle.n_pgrid * iz; //AWS20110202
	    
            for (ip = 0; ip < particle.n_pgrid  ; ++ip)                                              //AWS20110202
            {
	     Work_arr[ixyzp]        += alpha3_y_arr[ixyzp]        * particle_cr_density_arr[iixyzp]; //AWS20110202

                      ixyzp++;                                                                       //AWS20110202
	             iixyzp++;                                                                       //AWS20110202

	    }//ip               	    
	   }//iz
	  }//iy
 	 }//ix


	  // z direction -----------------------------------------------------------
             	  
#pragma omp parallel for schedule(runtime)  private(ix,iy,iz,ip,ixyzp,iixyzp) default(shared)       //AWS20110202
	   
	 for   (ix = 0; ix < particle.n_xgrid  ; ++ix)                                                //AWS20110202
         {
	  for  (iy = 0; iy < particle.n_ygrid  ; ++iy)                                               //AWS20110202
          {
  	   for (iz = 1; iz < particle.n_zgrid  ; ++iz)                                              //AWS20110202
           {
             ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid * iz;     //AWS20110202
	     iixyzp =  particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid *(iz-1);  //AWS20110202
	    
            for (ip = 0; ip < particle.n_pgrid  ; ++ip)                                              //AWS20110202
            {
	     Work_arr[ixyzp]        += alpha1_z_arr[ixyzp]        * particle_cr_density_arr[iixyzp]; //AWS20110202

                      ixyzp++;                                                                       //AWS20110202
	             iixyzp++;                                                                       //AWS20110202

	    }//ip               	    
	   }//iz
	  }//iy
 	 }//ix

#pragma omp parallel for schedule(runtime)  private(ix,iy,iz,ip,ixyzp,iixyzp) default(shared)       //AWS20110202
	   
	 for   (ix = 0; ix < particle.n_xgrid  ; ++ix)                                                //AWS20110202
         {
	  for  (iy = 0; iy < particle.n_ygrid  ; ++iy)                                               //AWS20110202
          {
  	   for (iz = 0; iz < particle.n_zgrid-1; ++iz)                                              //AWS20110202
           {
             ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid *  iz   ; //AWS20110202
	     iixyzp =  particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid * (iz+1); //AWS20110202
	    
            for (ip = 0; ip < particle.n_pgrid  ; ++ip)                                              //AWS20110202
            {
	     Work_arr[ixyzp]        += alpha3_z_arr[ixyzp]        * particle_cr_density_arr[iixyzp]; //AWS20110202

                      ixyzp++;                                                                       //AWS20110202
	             iixyzp++;                                                                       //AWS20110202

	    }//ip               	    
	   }//iz
	  }//iy
 	 }//ix


	  // p direction -----------------------------------------------------------
             	  
#pragma omp parallel for schedule(runtime)  private(ix,iy,iz,ip,ixyzp,iixyzp) default(shared)        //AWS20110202
	   
	 for   (ix = 0; ix < particle.n_xgrid  ; ++ix)                                               //AWS20110202
         {
	  for  (iy = 0; iy < particle.n_ygrid  ; ++iy)                                               //AWS20110202
          {
  	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)                                               //AWS20110202
           {
             ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid * iz     ;  //AWS20110202
	     iixyzp =  particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid * iz - 1 ;  //AWS20110202
	    
            for (ip = 1; ip < particle.n_pgrid  ; ++ip)                                              //AWS20110202
            {
	     Work_arr[ixyzp]        += alpha1_p_arr[ixyzp]        * particle_cr_density_arr[iixyzp]; //AWS20110202

                      ixyzp++;                                                                       //AWS20110202
	             iixyzp++;                                                                       //AWS20110202

	    }//ip               	    
	   }//iz
	  }//iy
 	 }//ix

#pragma omp parallel for schedule(runtime)  private(ix,iy,iz,ip,ixyzp,iixyzp) default(shared)       //AWS20110202
	   
	 for   (ix = 0; ix < particle.n_xgrid  ; ++ix)                                                //AWS20110202
         {
	  for  (iy = 0; iy < particle.n_ygrid  ; ++iy)                                               //AWS20110202
          {
  	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)                                              //AWS20110202
           {
             ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid *  iz     ; //AWS20110202
	     iixyzp =  particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid *  iz + 1 ; //AWS20110202
	    
            for (ip = 0; ip < particle.n_pgrid-1; ++ip)                                              //AWS20110202
            {
	     Work_arr[ixyzp]        += alpha3_p_arr[ixyzp]        * particle_cr_density_arr[iixyzp]; //AWS20110202

                      ixyzp++;                                                                       //AWS20110202
	             iixyzp++;                                                                       //AWS20110202

	    }//ip               	    
	   }//iz
	  }//iy
 	 }//ix



	 //---------------------------------------------------------------------

	   #pragma omp parallel for schedule(runtime)   default(shared)
	    for (ixyzp=0; ixyzp<nxyzp  ; ++ixyzp)                                                     //AWS20110202
	    {
	     particle_cr_density_arr[ixyzp]        +=     Work_arr[ixyzp];                            //AWS20110202
	    }



	  //-------------------------------------------------------------------------

	  // apply spatial boundary conditions
	  if( galdef.spatial_bound_conds==1) //AWS20131028
	  {
 
#pragma omp parallel for schedule(runtime) private(ix,iy,iz,ip) default(shared)
	  for (ip = 0; ip < particle.n_pgrid; ++ip) {
	    
	    for (ix = 0; ix < particle.n_xgrid; ++ix) {

	      for (iy = 0; iy < particle.n_ygrid; ++iy) {
	
		 iz = 0;  
	       	 ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid *  iz + ip    ;
		 particle_cr_density_arr[ixyzp] = 0;
		  
		 iz = particle.n_zgrid-1;
	     	 ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid *  iz + ip    ;
		 particle_cr_density_arr[ixyzp] = 0;

	      }
	     
	    }

     
	    for (ix = 0; ix < particle.n_xgrid; ++ix) {
	      
	      for (iz = 0; iz < particle.n_zgrid; ++iz) {
		    
		iy = 0;  
		ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid *  iz + ip    ;
		particle_cr_density_arr[ixyzp] = 0;

		
		
		iy = particle.n_ygrid-1;
		ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid *  iz + ip    ;
		particle_cr_density_arr[ixyzp] = 0;
	    

	      }	      

	    }
	     

	    for (iz = 0; iz < particle.n_zgrid; ++iz) {
	      
	      for (iy = 0; iy < particle.n_ygrid; ++iy) {
		
		ix = 0;  
		ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid *  iz + ip    ;
		particle_cr_density_arr[ixyzp] = 0;
	      		
		ix = particle.n_xgrid-1;
		ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid *  iz + ip    ;
		particle_cr_density_arr[ixyzp] = 0;
	      
		
	      }
	    
	    }

	  }  // ip  

	 } // spatial_bound_conds==1

	  //---------------------------------------------------------------------------------------------

// following needed only for diagnostics, so conditional for efficiency. 

	if(irept%galdef.timestep_diagnostics==0)
	  
	 for   (ix = 0; ix < particle.n_xgrid  ; ++ix)                                               //AWS20110202
         {
	  for  (iy = 0; iy < particle.n_ygrid  ; ++iy)                                               //AWS20110202
          {
  	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)                                              //AWS20110202
           {
            for (ip = 0; ip < particle.n_pgrid  ; ++ip)                                              //AWS20110202
            {
             ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid *  iz + ip    ; //AWS20110202
          
	     particle.cr_density.d3[ix][iy][iz].s[ip] =      particle_cr_density_arr[ixyzp];                    
	      
	    }//ip               	    
	   }//iz
	  }//iy
 	 }//ix


	 }//turbo==1                                                                                  //AWS20110202






	} // method == 2


	//--------------------------------------------------------------------------------------------------


	  if (0 == irept%galdef.timestep_print) {
	    
	    buf.str("");
	    buf<<"particle.cr_density after iteration "<<irept<<" with time step "<<dt/year2sec<<" yr";
	    INFO(buf.str());
	    particle.cr_density.print();
	    cout<<endl;
	  
	  }
	  
	  if (0 == irept%galdef.timestep_diagnostics) 
	  {
            converged =                                                //AWS20110114
	    propel_diagnostics(particle,
			       alpha1_x, alpha1_y, alpha1_z, alpha1_p,
			       alpha2_x, alpha2_y, alpha2_z, alpha2_p,
			       alpha3_x, alpha3_y, alpha3_z, alpha3_p,
			       total_source_function, 
			       dt);

            if(galdef.solution_convergence == 1 && converged == 1) //AWS20111020
	    {
	     buf.str("");
             buf<<" solution converged at irept="<<irept<<" exiting this timestep";
             INFO(buf.str());
	 
             break;// out of irept loop
            }
           }// if irept 

	}  //  irept

	if (2 == timestep_mode) 
	  timestep_mode2done = 1;
      
      }  //  dt>=end_timestep_sec
    
    }  //  timestep_mode
    

    //-------------------------------------------------------------------------------------------------------------

   // do this here since method 2 turbo does not do it in loop (for efficiency) except when diagnostics required

    if (galdef.solution_method == 21 && galdef.timestep_repeat2 > 0) //AWS20110202
	 for   (ix = 0; ix < particle.n_xgrid  ; ++ix)                                               //AWS20110202
         {
	  for  (iy = 0; iy < particle.n_ygrid  ; ++iy)                                               //AWS20110202
          {
  	   for (iz = 0; iz < particle.n_zgrid  ; ++iz)                                               //AWS20110202
           {
            for (ip = 0; ip < particle.n_pgrid  ; ++ip)                                              //AWS20110202
            {
             ixyzp =   particle.n_ygrid * particle.n_zgrid * particle.n_pgrid *  ix     + particle.n_zgrid * particle.n_pgrid *  iy    +  particle.n_pgrid *  iz + ip    ; //AWS20110202
          
	     particle.cr_density.d3[ix][iy][iz].s[ip] =      particle_cr_density_arr[ixyzp];                    
	      
	    }//ip               	    
	   }//iz
	  }//iy
 	 }//ix

   //-------------------------------------------------------------------------------------------------------------

    //Delete temporary arrays         //Gulli20070810
    alpha1_x.delete_array();
    alpha1_y.delete_array();
    alpha1_z.delete_array();
    alpha1_p.delete_array();
    alpha2_x.delete_array();
    alpha2_y.delete_array();
    alpha2_z.delete_array();
    alpha2_p.delete_array();
    alpha3_x.delete_array();
    alpha3_y.delete_array();
    alpha3_z.delete_array();
    alpha3_p.delete_array();
    
    Nx1_.delete_array();
    Ny1_.delete_array();
    Nz1_.delete_array();
    Np1_.delete_array();
    Nx2_.delete_array();
    Ny2_.delete_array();
    Nz2_.delete_array();
    Np2_.delete_array();
    Nx3_.delete_array();
    Ny3_.delete_array();
    Nz3_.delete_array();
    Np3_.delete_array();
    
    total_source_function.delete_array();
    Work.delete_array();


   if(galdef.solution_method==21)
   {
    delete[]      particle_cr_density_arr;                                                      //AWS20110204
    delete[]    total_source_function_arr;                                                      //AWS20110204
    delete[]                 alpha1_x_arr;                                                      //AWS20110204
    delete[]                 alpha1_y_arr;                                                      //AWS20110204
    delete[]                 alpha1_z_arr;                                                      //AWS20110204
    delete[]                 alpha1_p_arr;                                                      //AWS20110204
    delete[]                 alpha2_x_arr;                                                      //AWS20110204
    delete[]                 alpha2_y_arr;                                                      //AWS20110204
    delete[]                 alpha2_z_arr;                                                      //AWS20110204
    delete[]                 alpha2_p_arr;                                                      //AWS20110204
    delete[]                 alpha3_x_arr;                                                      //AWS20110204
    delete[]                 alpha3_y_arr;                                                      //AWS20110204
    delete[]                 alpha3_z_arr;                                                      //AWS20110204
    delete[]                 alpha3_p_arr;                                                      //AWS20110204
    delete[]                     Work_arr;                                                      //AWS20110204
   }

  } //3D case
  
  INFO("Exit");
  
  //cout<<"\n<<<<propel"<<endl;
  return 0;
}
