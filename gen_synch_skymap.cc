
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_synch_skymap.cc *                         galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// generate  synchrotron skymap 

#include <cassert>
#include <string>
#include <cstring>
#include "kappa_free_free.h"

using namespace std;//AWS20050624


#include "galprop_classes.h"
#include "galprop_internal.h"

#include <ErrorLogger.h>

double interpDistribution
 (Distribution &distribution,
  double x_min,double y_min,double z_min,
  double dx, double dy, double dz, 
  double x,double y,double z, int ip);                //AWS20110628

////////////////////////////////////////////////////////////////////////////////////

int Galprop::gen_synch_skymap() {

  int stat=0;
  INFO("Entry");
   
  double rad_to_deg=180./acos(-1.0); //AWS20110919

  if (galdef.skymap_format == 3 || 4 == galdef.skymap_format){
#pragma omp parallel for default(shared)
    for (int ip = 0; ip < galaxy.synchrotron_hp_skymap.Npix(); ++ip){
      SM::Coordinate co(galaxy.synchrotron_hp_skymap.pix2ang(ip));
      double l=co.l();
      double b=co.b();
      vector<double> synch, Q, U, free_free; //AWS20110905
      
      //     gen_synch_skymap_pixel(l, b,synch);
      gen_synch_skymap_pixel(l, b,synch,Q,U,free_free); //AWS20110905
      


      for(int inusynch=0; inusynch<galaxy.n_nu_synchgrid; inusynch++)
	{
	  galaxy.synchrotron_hp_skymap       [co][inusynch] =      synch [inusynch];
	  galaxy.synchrotron_Q_hp_skymap     [co][inusynch] =      Q     [inusynch]; //AWS20100709
	  galaxy.synchrotron_U_hp_skymap     [co][inusynch] =      U     [inusynch]; //AWS20100709
	  galaxy.synchrotron_P_hp_skymap     [co][inusynch] =      sqrt(Q[inusynch]*Q[inusynch] + U[inusynch]*U[inusynch]) ;//AWS20110328
	  galaxy.synchrotron_polang_hp_skymap[co][inusynch] = 0.5*atan2(U[inusynch],Q[inusynch]) * rad_to_deg ;//AWS20110919
	  galaxy.synchrotron_polfra_hp_skymap[co][inusynch] =      sqrt(Q[inusynch]*Q[inusynch] + U[inusynch]*U[inusynch])   /synch[inusynch];           //AWS20110922

	  galaxy.free_free_hp_skymap    [co][inusynch] = free_free[inusynch];//AWS20110905
	}
    }
    
    galaxy.synchrotron_hp_skymap       .setSpectra(&galaxy.nu_synch[0], galaxy.n_nu_synchgrid);
    galaxy.synchrotron_Q_hp_skymap     .setSpectra(&galaxy.nu_synch[0], galaxy.n_nu_synchgrid); //AWS20100709
    galaxy.synchrotron_U_hp_skymap     .setSpectra(&galaxy.nu_synch[0], galaxy.n_nu_synchgrid); //AWS20100709
    galaxy.synchrotron_P_hp_skymap     .setSpectra(&galaxy.nu_synch[0], galaxy.n_nu_synchgrid); //AWS20110328
    galaxy.synchrotron_polang_hp_skymap.setSpectra(&galaxy.nu_synch[0], galaxy.n_nu_synchgrid); //AWS20110919
    galaxy.synchrotron_polfra_hp_skymap.setSpectra(&galaxy.nu_synch[0], galaxy.n_nu_synchgrid); //AWS20110922

    galaxy.free_free_hp_skymap         .setSpectra(&galaxy.nu_synch[0], galaxy.n_nu_synchgrid); //AWS20110905
    
    if(galdef.verbose>=2)
      {
	cout<<" synchrotron skymap " <<endl;
	galaxy.synchrotron_hp_skymap.print(cout);
      }//galdef.verbose>=2
  }else{
#pragma omp parallel for schedule(dynamic) default(shared)
    for(int i_long=0; i_long<galaxy.n_long; i_long++)
      {
	for(int i_lat =0; i_lat <galaxy.n_lat; i_lat++)
	  {
	    double l=galaxy.long_min + i_long*galaxy.d_long;
	    double b=galaxy. lat_min + i_lat *galaxy.d_lat ;
	    vector<double> synch,Q,U,free_free;//AWS20110905
	    
	    //     gen_synch_skymap_pixel(l, b,synch);
	    gen_synch_skymap_pixel(l, b,synch,Q,U,free_free); //AWS20110905
	    
	    for(int inusynch=0; inusynch<galaxy.n_nu_synchgrid; inusynch++)
	      {
		galaxy.synchrotron_skymap       .d2[i_long][i_lat].s[inusynch] =      synch [inusynch];
		galaxy.synchrotron_Q_skymap     .d2[i_long][i_lat].s[inusynch] =      Q     [inusynch]; //AWS20100709
		galaxy.synchrotron_U_skymap     .d2[i_long][i_lat].s[inusynch] =      U     [inusynch]; //AWS20100709
		galaxy.synchrotron_P_skymap     .d2[i_long][i_lat].s[inusynch] =      sqrt(Q[inusynch]*Q[inusynch] + U[inusynch]*U[inusynch]) ;//AWS20110328
		galaxy.synchrotron_polang_skymap.d2[i_long][i_lat].s[inusynch] = 0.5*atan2(U[inusynch],Q[inusynch]) * rad_to_deg ;//AWS20110919
		galaxy.synchrotron_polfra_skymap.d2[i_long][i_lat].s[inusynch] =   sqrt(Q[inusynch]*Q[inusynch] + U[inusynch]*U[inusynch])   /synch[inusynch];           //AWS20110922

		galaxy.free_free_skymap    .d2[i_long][i_lat].s[inusynch] = free_free [inusynch];//AWS20110905
	      }
	  }
      }
    if(galdef.verbose>=2)
      {
	cout<<" synchrotron skymap " <<endl;
	galaxy.synchrotron_skymap.print();
      }//galdef.verbose>=2
  }
  INFO("Exit");
  return stat;
}
///////////////////////////////////////////////////////////////////////////////////////////
int Galprop::gen_synch_skymap_pixel(const double l, const double b, vector<double> &synch){
  double Ro= 8.5; // Galactocentric distance of Sun, kpc
  Ro= 8.3; // to avoid discontinuity in maps
  double dd0    =0.1; // max integration step in kpc at b=90 deg.
  double ddmax = 0.5; // max integration step allowed
  double dtr=acos(-1.)/180.;
  int ir,ix,iy,iz;
  //cout<<"l b ="<<l<<" "<<b<<endl;
  double sinb=sin(b*dtr);
  double cosb=cos(b*dtr);
  double sinl=sin(l*dtr);
  double cosl=cos(l*dtr);
  double d=0;
  double dd=dd0/(fabs(sinb)+1.e-6);
  if(dd>ddmax)dd=ddmax;

  dd = galdef.LoS_step ; //AWS20101103

  int complete=0;
  synch.resize(galaxy.n_nu_synchgrid, 0);
  while(complete==0)
    {
      d += dd;
      double RR=sqrt(Ro*Ro+pow(d*cosb,2)-2.0*Ro*d*cosb*cosl);
      double zz=d*sinb;
      double costheta=(Ro-d*cosb*cosl)/RR;
      if(costheta> 1.0)costheta= 1.0;
      if(costheta<-1.0)costheta=-1.0;
      
      if(gcr[0].n_spatial_dimensions==2)
	{
	  ir=(int)((RR-galaxy.r_min)/galaxy.dr + 0.5);//IMOS20060420
	  if(ir>galaxy.n_rgrid-1) { complete=1; ir=galaxy.n_rgrid-1; }
	  iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5);//IMOS20060420
	  if(iz<0               ) { complete=1; iz=0;                } 
	  if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }
	  // cout<<"d RR zz ir iz "<<d<<" "<<RR<<" "<<zz<<" "<<ir<<" "<<iz<<endl;
	} // particle.n_spatial_dimensions==2
      
      if(gcr[0].n_spatial_dimensions==3)
	{
	  int LRH=2; // 1=LH, 2=RH system                                                AWS20080311 AWS20101109 LH->RH
	  
	  double xx=Ro-d*cosb*cosl; // Sun on x axis at x=+Ro
	  double yy;
	  
	  if(LRH==1)                                                                 // AWS20080311  AWS20101109
	    yy=  +d*cosb*sinl; // Sun at y=0; +ve long=+ve y for Z=-X^Y (LH) system  // AWS20080311

          if(LRH==2)                                                                 // AWS20080311  AWS20101109
	    yy=  -d*cosb*sinl; // Sun at y=0; +ve long=-ve y for Z= X^Y (RH) system  // AWS20080311
	       
	       
               if(galdef.use_symmetry==1)
		 {
                  xx=fabs(xx);
                  yy=fabs(yy);
                  zz=fabs(zz);
               }
               ix=(int)((xx-galaxy.x_min)/galaxy.dx + 0.5);//IMOS20060420
               iy=(int)((yy-galaxy.y_min)/galaxy.dy + 0.5);//IMOS20060420
               iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5);//IMOS20060420

               if(ix<0               ) { complete=1; ix=0;                }
               if(iy<0               ) { complete=1; iy=0;                }  
               if(iz<0               ) { complete=1; iz=0;                } 
               if(ix>galaxy.n_xgrid-1) { complete=1; ix=galaxy.n_xgrid-1; }
               if(iy>galaxy.n_ygrid-1) { complete=1; iy=galaxy.n_ygrid-1; }
               if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }
//  cout<<"d RR xx yy zz ix iy iz "<<d<<" "<<RR<<" "<<xx<<" "<<yy<<" "<<zz<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
            } // particle.n_spatial_dimensions==3

            for(int inusynch=0; inusynch<galaxy.n_nu_synchgrid; inusynch++)
            {
              // float delta;   AWS20100810
	        double delta; //AWS20100810
               if(gcr[0].n_spatial_dimensions==2)
                  delta =dd*kpc2cm *galaxy.synchrotron_emiss.d2[ir][iz].    s[inusynch];

               if(gcr[0].n_spatial_dimensions==3)
		  delta =dd*kpc2cm *galaxy.synchrotron_emiss.d3[ix][iy][iz].s[inusynch];

               synch[inusynch] += delta;

//cout<<"l b ir iz  E_gamma bremss_ionized_emiss "<<l<<" "<<b<<" "<<ir<<" "<<iz<<" "<<" "
//<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.bremss_ionized_emiss.d2[ir][iz].s[iEgamma]<<endl;     
            }//inusynch
         }//complete==0
	return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////// Stokes version, AWS20100709
int Galprop::gen_synch_skymap_pixel(const double l, const double b, vector<double> &synch, vector<double> &Q,  vector<double> &U, vector<double>&free_free ) //AWS20110905
{
  double Ro= 8.5; // Galactocentric distance of Sun, kpc
  Ro= 8.3; // to avoid discontinuity in maps
  double dd0    =0.1; // max integration step in kpc at b=90 deg.
  double ddmax = 0.5; // max integration step allowed
  double dtr=acos(-1.)/180.;
  int ir,ix,iy,iz;
  //cout<<"l b ="<<l<<" "<<b<<endl;
  double sinb=sin(b*dtr);
  double cosb=cos(b*dtr);
  double sinl=sin(l*dtr);
  double cosl=cos(l*dtr);
  double d=0;
  double dd=dd0/(fabs(sinb)+1.e-6);
  if(dd>ddmax)dd=ddmax;

  dd = galdef.LoS_step ; //AWS20101103

  int integration_method = 2;//AWS20110629  1=original nearest gridpoint, 2=linear interpolation 
  if (galdef.verbose==-4000) integration_method=1;

  
  valarray<double> tau_free_free(galaxy.n_nu_synchgrid); //AWS20110630
                   tau_free_free = 0.0;                  //AWS20110630

  int complete=0;
  synch    .resize(galaxy.n_nu_synchgrid, 0);
  Q        .resize(galaxy.n_nu_synchgrid, 0);
  U        .resize(galaxy.n_nu_synchgrid, 0);
  free_free.resize(galaxy.n_nu_synchgrid, 0);            //AWS20110905       

  while(complete==0)
    {
      d += dd;
      double RR=sqrt(Ro*Ro+pow(d*cosb,2)-2.0*Ro*d*cosb*cosl);
      double zz=d*sinb;
      double costheta=(Ro-d*cosb*cosl)/RR;
      if(costheta> 1.0)costheta= 1.0;
      if(costheta<-1.0)costheta=-1.0;

      double xx,yy; //AWS20110628 moved from inside loop since required in interpolation
      
      if(gcr[0].n_spatial_dimensions==2)
	{
	  ir=(int)((RR-galaxy.r_min)/galaxy.dr + 0.5);//IMOS20060420
	  if(ir>galaxy.n_rgrid-1) { complete=1; ir=galaxy.n_rgrid-1; }
	  iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5);//IMOS20060420
	  if(iz<0               ) { complete=1; iz=0;                } 
	  if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }
	  // cout<<"d RR zz ir iz "<<d<<" "<<RR<<" "<<zz<<" "<<ir<<" "<<iz<<endl;
	} // particle.n_spatial_dimensions==2
      
      if(gcr[0].n_spatial_dimensions==3)
	{
	  int LRH=2; // 1=LH, 2=RH system                                                     AWS20080311 AWS20101109 LH->RH
	  
	   xx=Ro-d*cosb*cosl; // Sun on x axis at x=+Ro  AWS20110628
	                                    //double yy; AWS20110628
	  
	  if(LRH==1)                                                                       // AWS20080311  AWS20101109
	    yy=  +d*cosb*sinl; // Sun at y=0; +ve long=+ve y for Z=-X^Y (LH) system        // AWS20080311
          if(LRH==2)                                                                       // AWS20080311
	    yy=  -d*cosb*sinl; // Sun at y=0; +ve long=-ve y for Z= X^Y (RH) system        // AWS20080311  AWS20101109
	       
	       
               if(galdef.use_symmetry==1)
		 {
                  xx=fabs(xx);
                  yy=fabs(yy);
                  zz=fabs(zz);
               }
               ix=(int)((xx-galaxy.x_min)/galaxy.dx + 0.5);//IMOS20060420
               iy=(int)((yy-galaxy.y_min)/galaxy.dy + 0.5);//IMOS20060420
               iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5);//IMOS20060420

               if(ix<0               ) { complete=1; ix=0;                }
               if(iy<0               ) { complete=1; iy=0;                }  
               if(iz<0               ) { complete=1; iz=0;                } 
               if(ix>galaxy.n_xgrid-1) { complete=1; ix=galaxy.n_xgrid-1; }
               if(iy>galaxy.n_ygrid-1) { complete=1; iy=galaxy.n_ygrid-1; }
               if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }
//  cout<<"d RR xx yy zz ix iy iz "<<d<<" "<<RR<<" "<<xx<<" "<<yy<<" "<<zz<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
            } // particle.n_spatial_dimensions==3



            for(int inusynch=0; inusynch<galaxy.n_nu_synchgrid; inusynch++)
            {
	    
	      double delta(0),deltaQ(0),deltaU(0),delta_free_free(0); // AWS20110905

               if(gcr[0].n_spatial_dimensions==2)
	       {
                  delta =dd*kpc2cm *galaxy.synchrotron_emiss  .d2[ir][iz].    s[inusynch];
                  deltaQ=dd*kpc2cm *galaxy.synchrotron_Q_emiss.d2[ir][iz].    s[inusynch];
                  deltaU=dd*kpc2cm *galaxy.synchrotron_U_emiss.d2[ir][iz].    s[inusynch];
	       }

               if(gcr[0].n_spatial_dimensions==3)
	       {

                  if(integration_method==1)//AWS20110629
		  {

		   delta =dd*kpc2cm *galaxy.synchrotron_emiss  .d3[ix][iy][iz].s[inusynch];
		   deltaQ=dd*kpc2cm *galaxy.synchrotron_Q_emiss.d3[ix][iy][iz].s[inusynch];
		   deltaU=dd*kpc2cm *galaxy.synchrotron_U_emiss.d3[ix][iy][iz].s[inusynch];

                  }//integration_method==1


                  if(integration_method==2)//AWS20110629
		  {

		   delta=dd*kpc2cm*interpDistribution(galaxy.synchrotron_emiss,galaxy. x_min,galaxy .y_min,galaxy .z_min,
						                               galaxy.dx,    galaxy.dy,    galaxy.dz,
                                                                                      xx,           yy,           zz,
                                                      inusynch);

		   deltaQ=dd*kpc2cm*interpDistribution(galaxy.synchrotron_Q_emiss,galaxy. x_min,galaxy .y_min,galaxy .z_min,
						                                  galaxy.dx,    galaxy.dy,    galaxy.dz,
                                                                                         xx,           yy,           zz,
                                                      inusynch);

		   deltaU=dd*kpc2cm*interpDistribution(galaxy.synchrotron_U_emiss,galaxy. x_min,galaxy .y_min,galaxy .z_min,
						                                  galaxy.dx,    galaxy.dy,    galaxy.dz,
                                                                                         xx,           yy,           zz,
                                                      inusynch);


                  }//integration_method==2
	       } //3D



               //cout<<"l b ir iz  E_gamma bremss_ionized_emiss "<<l<<" "<<b<<" "<<ir<<" "<<iz<<" "<<" "
               //<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.bremss_ionized_emiss.d2[ir][iz].s[iEgamma]<<endl; 
    
               //double Ne=1,Te=7000; just for testing
	       //kappa_free_free( galaxy.nu_synch[inusynch], Ne,  Te, 0,1);
	       //kappa_free_free( galaxy.nu_synch[inusynch], xx,yy,zz, Te, 0,1);

               if(galdef.free_free_absorption>=1)
	       {
		 double Te=galdef.HII_Te;
                 double clumping_factor=galdef.HII_clumping_factor;//AWS20110704
                 double emiss_free_free;

                int debug=0;

                if(gcr[0].n_spatial_dimensions==2)
		  tau_free_free[inusynch] += 
                   dd*kpc2cm*kappa_free_free_2D( galaxy.nu_synch[inusynch], RR,   zz, Te, clumping_factor,emiss_free_free,0,debug);//AWS20110704

                if(gcr[0].n_spatial_dimensions==3)
		  tau_free_free[inusynch] += 
                   dd*kpc2cm*kappa_free_free_3D( galaxy.nu_synch[inusynch], xx,yy,zz, Te, clumping_factor,emiss_free_free,0,debug);//AWS20110704

                delta_free_free = dd*kpc2cm*emiss_free_free; //AWS20110905

		// replace total  synchrotron with free-free 
                if(galdef.free_free_absorption==2){ delta=dd*kpc2cm*emiss_free_free;  }

                double absorption=exp(-tau_free_free[inusynch]);


	        if(debug==1)
	        {
                 cout<<"nu="<< galaxy.nu_synch[inusynch]<<endl;
                 cout<< " tau_free_free="<<tau_free_free[inusynch]<<endl;
                 cout<< " absorption ="<<absorption<<endl;
	        } 

		delta           *=absorption;
                deltaQ          *=absorption;
                deltaU          *=absorption;
		delta_free_free *=absorption;                //AWS20110905

	       }//absorption

               synch    [inusynch] += delta;
               Q        [inusynch] += deltaQ;
               U        [inusynch] += deltaU;
               free_free[inusynch] += delta_free_free;       //AWS20110905

            }//inusynch

         }//complete==0

	return 0;
}
///////////////////////////////////////////////////////////////AWS20110628

double interpDistribution
 (Distribution &distribution,
  double  x_min,double  y_min,double  z_min,
  double dx,    double dy,    double dz, 
  double  x,    double  y,    double  z,
  int ip)
{
  // 3D linear interpolation at (x,y,z) for one spectral point
  // successive interpolations, first in x, then y, then z
  // it will eventually be turned into a spline version

          if(distribution.n_spatial_dimensions!=3) return -1.0;

          int debug=0;
 
          int     ix=(int)((x-x_min)/dx + 0.0001); // grid point  below interpolation point
          int     iy=(int)((y-y_min)/dy + 0.0001);
          int     iz=(int)((z-z_min)/dz + 0.0001);

	  // stay within grid for all cases
          if(ix<0)                     ix=0;
          if(iy<0)                     iy=0;
          if(iz<0)                     iz=0;
          if(ix>distribution.n_xgrid-2)ix=distribution.n_xgrid-2;
          if(iy>distribution.n_ygrid-2)iy=distribution.n_ygrid-2;
          if(iz>distribution.n_zgrid-2)iz=distribution.n_zgrid-2;

          double x1=x_min + ix*dx;                 // coordinates of lower grid point
          double y1=y_min + iy*dy;
          double z1=z_min + iz*dz;

          double v1=distribution.d3[ix  ][iy]  [iz  ].s[ip];
          double v2=distribution.d3[ix+1][iy]  [iz  ].s[ip];
          double v3=distribution.d3[ix  ][iy+1][iz  ].s[ip];
          double v4=distribution.d3[ix+1][iy+1][iz  ].s[ip];
          double v5=distribution.d3[ix  ][iy  ][iz+1].s[ip];
          double v6=distribution.d3[ix+1][iy  ][iz+1].s[ip];
          double v7=distribution.d3[ix  ][iy+1][iz+1].s[ip];
          double v8=distribution.d3[ix+1][iy+1][iz+1].s[ip];



          double v12= v1  +  (v2-v1)  *(x-x1)/dx; // iy   iz   at x
          double v34= v3  +  (v4-v3)  *(x-x1)/dx; // iy+1 iz   at x
          double v56= v5  +  (v6-v5)  *(x-x1)/dx; // iy   iz+1 at x
          double v78= v7  +  (v8-v7)  *(x-x1)/dx; // iy+1 iz+1 at x

          double w1 = v12 +  (v34-v12)*(y-y1)/dy; //      iz   at x,y
          double w2 = v56 +  (v78-v56)*(y-y1)/dy; //      iz+1 at x,y


	  double r  = w1  +  (w2-w1)  *(z-z1)/dz;//   at x,y,z

          if(debug==1)
	  {
	      cout    <<endl<<"interpDistribution:"<<endl;
	      cout    << " x=" << x<<  " y="<<  y<< " z="<<  z <<" ip="<<ip<<endl;
	      cout    <<" ix=" <<ix<< " iy="<< iy<<" iz="<< iz<<endl;
	      cout    << " x1="<< x1<< " y1="<< y1<<" z1="<< z1<<endl;
	      cout    << " v1="<< v1<< " v2="<< v2<<" v3="<< v3 << " v4="<< v4<<endl;
	      cout    << " v5="<< v5<< " v6="<< v6<<" v7="<< v7 << " v8="<< v8<<endl;
	      cout    << " v12="<< v12<< " v34="<< v34<<" v56="<< v56 << " v78="<< v78<<endl;
	      cout    << " w1="<< w1<< " w2="<<w2<<endl;
              cout    << " r="<< r<<endl;
	      cout    <<endl;
	  }

          return r;
}
