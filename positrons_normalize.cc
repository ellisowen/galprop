
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * positrons_normalize.cc *                      galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// based on electrons_normalize.cc AWS20131023

using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"
#include <cstring>

#include "ErrorLogger.h"
#include <sstream>

int Galprop::positrons_normalize()
{
   INFO("Entry");

// identify the primary positrons 
  int ielectrons=-1;               // variable name kept from electrons_normalize.cc
  for(int i=0;i<n_species;i++) 
    if(strcmp(gcr[i].name,"primary_positrons")==0) ielectrons=i;
  if(ielectrons==-1){WARNING("primary positrons not found!"); return 1;}
  ostringstream buf;
  buf<<"  primary positrons found as species #"<<ielectrons;
  INFO(buf.str());

  double  r0=8.5; // solar Galactocentric radius, kpc
// deriving the normalization grid point  
  int ip=(int)(log(galdef.electron_norm_Ekin/galdef.Ekin_min)/log(galdef.Ekin_factor) + 1.0e-6);
  int iz=(int)((0.-galdef.z_min)                             /galdef.dz               + 1.0e-6); // z=0, Galactic plane                   
// check normalization point in Ekin  AWS20131107 NB galaxy.n_pgrid not defined!
  if(ip<0 || ip+1>gcr[ielectrons].n_pgrid-1   )
  {
    buf.str("");
    buf<<"ERROR: positrons normalization Ekin: ip="<<ip<<": "<< galdef.positron_norm_Ekin << " outside range of grid "
       <<  gcr[ielectrons].Ekin[0] <<" - "<<  gcr[ielectrons].Ekin[gcr[ielectrons].n_pgrid-1]
       <<" gcr[ielectron].n_pgrid="<<gcr[ielectrons].n_pgrid
       << " ; no normalization can be done, exiting!";
    INFO(buf.str());
    exit(1);
  }    

  double v1,v2,v3,v4,v5,v6;
  if(galdef.n_spatial_dimensions==2)
    {
      int ir=(int)((r0-galdef.r_min)/galdef.dr + 1.0e-6);
      buf.str("");
      buf<<"Grid point for normalization: ir r[ir] iz z[iz] ip Ekin[ip]:Ekin[ip+1] "<<ir<<" " <<gcr[ielectrons].r[ir]<<" " <<iz <<" "<<gcr[ielectrons].z[iz]
         <<" "<<ip<<" "<< gcr[ielectrons].Ekin[ip]<<":"<< gcr[ielectrons].Ekin[ip+1]
          <<" normalization Ekin="<<galdef.positron_norm_Ekin ;   //AWS20131107

      INFO(buf.str());
      
      v1=gcr[ielectrons].cr_density.d2[ir  ][iz].s[ip];
      v2=gcr[ielectrons].cr_density.d2[ir+1][iz].s[ip];
      v3=gcr[ielectrons].cr_density.d2[ir  ][iz].s[ip+1];
      v4=gcr[ielectrons].cr_density.d2[ir+1][iz].s[ip+1];
      v5=v1+(r0-gcr[ielectrons].r[ir])/galdef.dr*(v2-v1); // r0 ip
      v6=v3+(r0-gcr[ielectrons].r[ir])/galdef.dr*(v4-v3); // r0 ip+1
    }//n_spatial_dimensions==2
  
  if(galdef.n_spatial_dimensions==3)
    {
      int ix=(int)((r0-galdef.x_min)/galdef.dx + 1.0e-6);
      int iy=(int)((0.-galdef.y_min)/galdef.dy + 1.0e-6); 
      buf.str("");
       buf<<"Grid point for normalization: ix x[ix] iy y[iy] iz z[iz] ip Ekin[ip]:Ekin[ip+1] "
          <<ix<<" " <<gcr[ielectrons].x[ix]
          <<" "<<iy<<" "<< gcr[ielectrons].y[iy]<<" " <<iz <<" "<<gcr[ielectrons].z[iz]
          <<" "<<ip<<" "<< gcr[ielectrons].Ekin[ip]<<":"<< gcr[ielectrons].Ekin[ip+1]
          <<" normalization Ekin="<<galdef.positron_norm_Ekin ;   //AWS20131107
      INFO(buf.str());
      
      v1=gcr[ielectrons].cr_density.d3[ix  ][iy][iz].s[ip];   //AWS20001121
      v2=gcr[ielectrons].cr_density.d3[ix+1][iy][iz].s[ip];   //AWS20001121
      v3=gcr[ielectrons].cr_density.d3[ix  ][iy][iz].s[ip+1]; //AWS20001121
      v4=gcr[ielectrons].cr_density.d3[ix+1][iy][iz].s[ip+1]; //AWS20001121
      v5=v1+(r0-gcr[ielectrons].x[ix])/galdef.dx*(v2-v1); // r0 ip
      v6=v3+(r0-gcr[ielectrons].x[ix])/galdef.dx*(v4-v3); // r0 ip+1
      
    }//n_spatial_dimensions==3
  
  double vnorm=exp( log(v5)+log(galdef.positron_norm_Ekin/gcr[ielectrons].Ekin[ip])/log(galdef.Ekin_factor)*log(v6/v5) );
  
  buf.str("");
  buf<<"v1 v2 v3 v4 v5 v6 vnorm  "<<v1<<" " <<v2<<" " <<v3 <<" "<<v4 <<" "<<v5<<" "<< v6<<" "<<vnorm;
  INFO(buf.str());

 // normalize positrons
  galdef.electron_source_normalization *= galdef.positron_norm_flux/vnorm; // IMOS20031016
  gcr[ielectrons].cr_density           *=(galdef.positron_norm_flux/vnorm);
  gcr[ielectrons].normalization_factor  =(galdef.positron_norm_flux/vnorm); //AWS20010121
  
  if(galdef.verbose>=2){
    INFO("primary positrons after normalization:");
    gcr[ielectrons].cr_density.print();
  }  

  INFO("Exit");
return 0;
}
