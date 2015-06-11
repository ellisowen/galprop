#include<iostream>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

void B_field_3D_model
(const std::string &name, const std::vector<double> &parameters,
 double x, double y, double z,
 int options,
 double &Breg,  double &Bregx, double &Bregy, double &Bregz,
 double &Bran,  double &Branx, double &Brany, double &Branz,
 int debug=0 )
{

  // General purpose interface for Galactic magnetic field models
  // including regular and random components

  // input  arguments

  // name: name of model
  // parameters: array of model parameters for model with this name
  // options:    for addional control
  // x,y,z  Galactic coordinates of point, kpc  (Sun at x=Rsun, y=z=0)


  // output arguments
  // magnetic field strengths in Gauss
  // Breg   regular field total
  // Bregx  regular field x-component
  // Bregy  regular field y-component
  // Bregz  regular field z-component

  // Bran   random  field total
  // Branx  random  field x-component (for a random direction)
  // Brany  random  field y-component (for a random direction) 
  // Branz  random  field z-component (for a random direction)

  // local variables
  double theta,phi; // spherical angles
  double dtr=acos(-1.)/180.; // degrees to radians
  double pi =acos(-1.);

  if(debug==1)cout<<"B_field_3D_model >>"<<endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  if(name == "test")
  {
    // uniform regular field strength and direction, and uniform random field
    Breg=6e-6;
    theta=89*dtr;
    phi= 100*dtr;
    Bregx=Breg*sin(theta)*cos(phi);
    Bregy=Breg*sin(theta)*sin(phi);
    Bregz=Breg*cos(theta);         

    Bran=5e-6;
    theta=30*dtr;
    phi= 200*dtr;
    Branx=Bran*sin(theta)*cos(phi);
    Brany=Bran*sin(theta)*sin(phi);
    Branz=Bran*cos(theta);  
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(name == "circular")
  {
    // uniform regular field strength parallel to plane, circular geometry, and uniform random field
    Breg=6e-6;
    theta=90.*dtr;  // zenith angle for in-plane field
    phi= atan2(y,x)+ pi/2; // this makes field direction circular
    Bregx=Breg*sin(theta)*cos(phi);
    Bregy=Breg*sin(theta)*sin(phi);
    Bregz=Breg*cos(theta);         

    Bran=5e-6;
    //   theta=30*dtr;  in future should assign a random direction
    //    phi= 200*dtr;
    Branx=Bran*sin(theta)*cos(phi);
    Brany=Bran*sin(theta)*sin(phi);
    Branz=Bran*cos(theta);  
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(name == "circular2")
  {
    // uniform regular field strength parallel to plane, circular geometry, and uniform random field
    // parameterized as below
    Breg=parameters[0];
    theta=90.*dtr;  // zenith angle for in-plane field
    phi= atan2(y,x)+ pi/2; // this makes field direction circular
    Bregx=Breg*sin(theta)*cos(phi);
    Bregy=Breg*sin(theta)*sin(phi);
    Bregz=Breg*cos(theta);         

    Bran=parameters[1];
    //   theta=30*dtr;  in future should assign a random direction
    //    phi= 200*dtr;
    Branx=Bran*sin(theta)*cos(phi);
    Brany=Bran*sin(theta)*sin(phi);
    Branz=Bran*cos(theta);  
  }



 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(name == "spiral")
  {
    // uniform regular field strength parallel to plane, spiral geometry, and uniform random field
    // parameterized as below
    Breg=parameters[0];
    double pitch_angle=parameters[1]*dtr; 
    theta=90.*dtr;  // zenith angle for in-plane field
    phi= atan2(y,x)+ pi/2 + pitch_angle; // pitch angle relative to circle
    Bregx=Breg*sin(theta)*cos(phi);
    Bregy=Breg*sin(theta)*sin(phi);
    Bregz=Breg*cos(theta);         

    Bran=parameters[2];
    //   theta=30*dtr;  in future should assign a random direction
    //    phi= 200*dtr;
    Branx=Bran*sin(theta)*cos(phi);
    Brany=Bran*sin(theta)*sin(phi);
    Branz=Bran*cos(theta);  
  }






  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  if(name == "galprop_original")
  {
    // same model as in galprop B_field_model.cc 
    // Bo rscale zscale encoded in 9-digit number: BBBrrrzzz in units of 0.1
    // e.g. 123456789 : Bo= 12.3E-10 Tesla  rscale=45.6 kpc zscale=78.9 kpc



    double Bo,rscale,zscale;
    double ro=8.5;
    double b_field;

    

    double r=sqrt(x*x+y*y);


    /* this was original test case, no longer used
    int model=int(parameters[0]+.001);// using this to represent galplotdef.B_field_model; parameters are double so convert to integer for galprop B_field_model.cc case
                                      // NB avoid leading zeros when defining B_field_model using integer (it is interpreted as hex!)
    // following 7-line  code segment is copied from  galprop B_field_model.cc, gives total B random in Tesla
   if (model > 1000)
   {
      Bo=           (model/1000000)                 * 0.1 *1.0e-10;
      rscale=(model-(model/1000000)*1000000 )/1000  * 0.1         ;
      zscale=(model%1000)                           * 0.1         ;
      b_field=Bo *exp(-(r-ro)/rscale) * exp(-fabs(z)/zscale);

   if(debug==1)
   cout<<"galprop original model="<<model<<" Bo="<<Bo<<" rscale="<<rscale<<" zscale="<<zscale
    <<" (x, y, z) = ("<<x<<", "<<y<<", "<<z<<")" << " r=" <<r 
    <<" b_field="<<b_field<<endl;

   }
    */

   // now take parameters as specified for 3D model in galdef file           AWS20080314
   Bo   = parameters[0]; // Gauss
   rscale=parameters[1]; // kpc
   zscale=parameters[2]; // kpc
   b_field=Bo *exp(-(r-ro)/rscale) * exp(-fabs(z)/zscale);

 

   // the regular field is zero
   Breg =0.0;
   Bregx=0.0;
   Bregy=0.0;
   Bregz=0.0;

   Bran=b_field;
   Branx=Bran;  // put all field in x-direction (has no significance)
   Brany=0.0;
   Branz=0.0;

   theta=0.0;
   phi  =0.0;

  } //galprop_original

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

    // WMAP regular field strength and direction (x,y,z); right handed system; Sun in x; B in Page et al. 2007 adapted to this coordinate system
    // from Elena Orlando: but she says this formulation is invalid: keep for reference
if(name == "wmap_page")
  {
  
       
  
   double chi0,chi,psi0,psi,psi1;
   double ro=8;

    Breg=parameters[0];
    phi=atan2(y,x);
    chi0=parameters[1]*dtr;//parameter[1]=25 in Page et al 2007; z dependence
    chi=chi0*tanh(z);
    psi0=parameters[2]*dtr;//perameter[2]=35; opening angle of spiral arms
    psi1=parameters[3]*dtr;//parameter[3]=0.9 radial dependence of opening angle
    psi=psi0+psi1*log(pow(x*x+y*y,0.5)/ro);
    Bregx=Breg*cos(chi)*sin(psi+phi);
    Bregy=-Breg*cos(chi)*cos(psi+phi);
    Bregz=Breg*sin(chi);  

   Bran=parameters[4]; 
   Branx=Bran;  // put all field in x-direction (has no significance)
   Brany=0.0;
   Branz=0.0;

   theta=0.0;


}

/////////////////////////////////////////////////////////////////////////////////////////////

// Han, J.L. 1994 A&A 288,759-772 formulation = WMAP Miville-Deschenes et al.  2008 arXiv:0802.3345
// select the parameters. bi-symmetric spiral; Sun on x-axis  , theta=clockwise 
// in the direction of positive y; reference system of galprop.
// A constant random component can be added via parameters[4]
// from Elena Orlando

if(name == "han")
 {

 
   //    Breg=parameters[0];           //=1.8 const in Han; =3 const in WMAP M.        //AWS20080331 Breg is output parameter!

    double B0=parameters[0];           //=1.8 const in Han; =3 const in WMAP M.        //AWS20080331
    theta=atan2(y,x);             // NB theta defined differently from above 
    double chi0=parameters[1]*dtr;// not defined in Han; parameter[1]=8 in WMAP Miville
    double chi=chi0*tanh(z);      //  z dependence, not defined in Han;
    double p=parameters[2]*dtr;   // pitch angle =-8.2 in Han;-8.5 in WMAP M.
    double psi=1./(tan(p));       // radial dependence of opening angle
    double ro=parameters[3];      // =11.9 in Han, =11 in WMAP M.     
    
    //    Bregx=Breg*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*sin(p-theta)*cos(chi);    //AWS20080331
    //    Bregy=Breg*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*cos(p-theta)*cos(chi);    //AWS20080331
    Bregx=B0*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*sin(p-theta)*cos(chi);            //AWS20080331
    Bregy=B0*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*cos(p-theta)*cos(chi);            //AWS20080331

    Bregz=0.0; // Han (and Miville-Deschenes ??) have zero z-component
    Bregz=parameters[5];                                                               //AWS20080711

    //Bregz=Breg*sin(chi); in case we need a z-component in future

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy +  Bregz*Bregz);                               //AWS20080331

    // constant random component can be added 
    Bran=parameters[4]; 
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

   
}
/////////////////////////////////////////////////////////////////////////////////////////////
// same as model "han" but random field is scaled to total regular at each point using parameters[4]
// contributed by AWS

  if(name == "han_scaled_ran") //AWS20080319
 {



    double B0=parameters[0];           //=1.8 const in Han; =3 const in WMAP M.        //AWS20080331
    theta=atan2(y,x);             // NB theta defined differently from above 
    double chi0=parameters[1]*dtr;// not defined in Han; parameter[1]=8 in WMAP Miville
    double chi=chi0*tanh(z/parameters[9]);      //  z dependence, not defined in Han;  //EO20081120
    //   double chi=tanh(z/parameters[9]);      //  z dependence, not defined in Han; error in Miville formula ??             //EO20081120
    double p=parameters[2]*dtr;   // pitch angle =-8.2 in Han;-8.5 in WMAP M.
    double psi=1./(tan(p));       // radial dependence of opening angle
    double ro=parameters[3];      // =11.9 in Han, =11 in WMAP M.     
    
  
    Bregx=B0*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*sin(p-theta)*cos(chi);            //AWS20080331
    Bregy=B0*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*cos(p-theta)*cos(chi);            //AWS20080331

    Bregz=0.0; // Han (and Miville-Deschenes ??) have zero z-component
    Bregz=parameters[5];                                                               //AWS20080711

    //Bregz=Breg*sin(chi); in case we need a z-component in future

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy +  Bregz*Bregz);                               //AWS20080331

    // random component proportional to total regular
    //    Bran=parameters[4]*sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);            //AWS20080331

    Bran=parameters[4]*Breg;                                                           //AWS20080331
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

   
}

/////////////////////////////////////////////////////////////////////////////////////////////
// same as model "han" but random field constant but with z dependence


  if(name == "han_ran") //EO20081127
 {



    double B0=parameters[0];           
    theta=atan2(y,x);             // NB theta defined differently from above 
    double chi0=parameters[1]*dtr;// not defined in Han; parameter[1]=8 in WMAP Miville
    double chi=chi0*tanh(z/parameters[9]);      //  z dependence, not defined in Han;  //EO20081120
    //   double chi=tanh(z/parameters[9]);      //  z dependence, not defined in Han; error in Miville formula ??             //EO20081120
    double p=parameters[2]*dtr;   // pitch angle =-8.2 in Han;-8.5 in WMAP M.
    double psi=1./(tan(p));       // radial dependence of opening angle
    double ro=parameters[3];      // =11.9 in Han, =11 in WMAP M.     
    
  
    Bregx=B0*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*sin(p-theta)*cos(chi);            //AWS20080331
    Bregy=B0*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*cos(p-theta)*cos(chi);            //AWS20080331

    Bregz=0.0; // Han (and Miville-Deschenes ??) have zero z-component
    Bregz=parameters[5];                                                              

    //Bregz=Breg*sin(chi); in case we need a z-component in future

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy +  Bregz*Bregz);                            

  

    Bran=parameters[4]*cos(chi);                                                         
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

   
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Han with exponential z dependence, new model
// from Elena Orlando

 if(name == "Han_exp_scaled")           
 {
    double B0=parameters[0];
                                                         
   double  ro=8.5;
   theta=atan2(y,x);
   double p=parameters[2]*dtr;//=-8 pitch angle
   double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   //   double z0=1.5;                                            
   double z0=parameters[6];                         
   //   double fz=(z/fabs(z))*exp(-fabs(z/z0));   
   double fz=exp(-fabs(z/z0));//                                    
   if(z<0.) fz=-fz; // to avoid problem at z=0
   double  psi=1./(tan(p));//radial dependence of opening angle   
   double r=pow(x*x+y*y,0.5);


    Bregx=B0  *cos(theta-psi*log(    r     /r0))*sin(p-theta)*fz;                        
    Bregy=B0  *cos(theta-psi*log(    r     /r0))*cos(p-theta)*fz;
    Bregz=0.;
   
    Bregz=parameters[5];                                                            


    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                                 

    // random component can be added for testing
    Bran=parameters[4]*Breg;
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }

//////////////////////////////////////////////////////////////////////////////////////////////////
// Han with exponential z dependence, new model   EO20081112
// from Elena Orlando

 if(name == "Han_exp")
 {
    double B0=parameters[0];

   double  ro=8.5;
   theta=atan2(y,x);
   double p=parameters[2]*dtr;//=-8 pitch angle
   double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   //   double z0=1.5;
   double z0=parameters[6];
   //   double fz=(z/fabs(z))*exp(-fabs(z/z0));
   double fz=exp(-fabs(z/z0));//
   if(z<0.) fz=-fz; // to avoid problem at z=0
   double  psi=1./(tan(p));//radial dependence of opening angle
   double r=pow(x*x+y*y,0.5);

    Bregx=B0  *cos(theta-psi*log(    r     /r0))*sin(p-theta)*fz;
    Bregy=B0  *cos(theta-psi*log(    r     /r0))*cos(p-theta)*fz;
    Bregz=0.;

    Bregz=parameters[5];

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);

    // random component can be added for testing
    Bran=parameters[4];
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }
/////////////////////////////////////////////////////////////////////////////////////////////////

// Han with exponential z dependence, new model   EO20081112
// from Elena Orlando

 if(name == "Han_ran_exp")
 {
    double B0=parameters[0];

   double  ro=8.5;
   theta=atan2(y,x);
   double p=parameters[2]*dtr;//=-8 pitch angle
   double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   //   double z0=1.5;
   double z0=parameters[6];
   //   double fz=(z/fabs(z))*exp(-fabs(z/z0));
   double fz=exp(-fabs(z/z0));//
   if(z<0.) fz=-fz; // to avoid problem at z=0
   double  psi=1./(tan(p));//radial dependence of opening angle
   double r=pow(x*x+y*y,0.5);

    Bregx=B0  *cos(theta-psi*log(    r     /r0))*sin(p-theta)*fz;
    Bregy=B0  *cos(theta-psi*log(    r     /r0))*cos(p-theta)*fz;
    Bregz=0.;

    Bregz=parameters[5];

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);

    // random component can be added for testing
    Bran=parameters[4];
    Branx=Bran*fz;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }


//////////////////////////////////////////////////////////////////////////////////////////////////

// Tinyakov & Tkachev 2002,Astroparticle Physics 18, 165  BSS-A model, radial dependence of Bo, no Bz
// from Elena Orlando

 if(name == "tinyakov")            //AWS20080311
 {
   double B0;                                                                              //AWS20080331
   double  ro=8.5;
   theta=atan2(y,x);
   double p=parameters[2]*dtr;//=-8 pitch angle
   double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   //   double z0=1.5;                                               AWS20081107
   double z0=parameters[6];                          //              AWS20081107
   //   double fz=(z/fabs(z))*exp(-fabs(z/z0));      //              AWS20080609
   double fz=exp(-fabs(z/z0));//                                     AWS20080609
   if(z<0.) fz=-fz; // to avoid problem at z=0
   double  psi=1./(tan(p));//radial dependence of opening angle    //AWS20080609
   double r=pow(x*x+y*y,0.5);

   double rbreak = parameters[8];                                  //EO20081119

    if(r<=rbreak) B0  =parameters[0]*(ro/rbreak);//parameters[0]=1.4
    if(r> rbreak) B0  =parameters[0]*(ro/r);

    Bregx=B0  *cos(theta-psi*log(    r     /r0))*sin(p-theta)*fz;                          //AWS20080311 20081107
    Bregy=B0  *cos(theta-psi*log(    r     /r0))*cos(p-theta)*fz; //correction AWS20080320   AWS20080311 20081107
    Bregz=0.;
    Bregz=parameters[5];                                                               //AWS20080711

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy +  Bregz*Bregz);                                   //AWS20080331

    // random component can be added for testing
    Bran=parameters[4]; 
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }







//////////////////////////////////////////////////////////////////////////////////////////////////

// Tinyakov & Tkachev 2002,Astroparticle Physics 18, 165  BSS-A model, radial dependence of Bo, no Bz
// from Elena Orlando
// same as model "tinyakov" but random field is scaled to total regular at each point using parameters[4]

 if(name == "tinyakov_scaled_ran")            //AWS20080311
 {
   double B0;                                                                              //AWS20080331
   double  ro=8.5;
   theta=atan2(y,x);
   double p=parameters[2]*dtr;//=-8 pitch angle
   double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   //   double z0=1.5;                                               AWS20081107
   double z0=parameters[6];                          //              AWS20081107
   //   double fz=(z/fabs(z))*exp(-fabs(z/z0));      //              AWS20080609
   double fz=exp(-fabs(z/z0));//                                     AWS20080609
   if(z<0.) fz=-fz; // to avoid problem at z=0
   double  psi=1./(tan(p));//radial dependence of opening angle    //AWS20080609
   double r=pow(x*x+y*y,0.5);

   double rbreak = parameters[8];                                  //EO20081119

    if(r<=rbreak) B0  =parameters[0]*(ro/rbreak);//parameters[0]=1.4
    if(r> rbreak) B0  =parameters[0]*(ro/r);

    Bregx=B0  *cos(theta-psi*log(    r     /r0))*sin(p-theta)*fz;                          //AWS20080311  20081107
    Bregy=B0  *cos(theta-psi*log(    r     /r0))*cos(p-theta)*fz; //correction AWS20080320   AWS20080311  20081107
    Bregz=0.;
    Bregz=parameters[5];                                                               //AWS20080711

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy +  Bregz*Bregz);                                   //AWS20080331

  


    Bran=parameters[4]*Breg;                                                           //AWS20080331
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;



 }

 //////////////////////////////////////////////////////////////////////////////////////
// Tinyakov & Tkachev 2002,Astroparticle Physics 18, 165  BSS-A model,
// radial dependence of Bo
// from Elena Orlando

 if(name == "tinyakov_exp_ran")            //EO20081118
 {
   double B0;                                                                             

   double  ro=8.5;
   theta=atan2(y,x);
   double p=parameters[2]*dtr;                      //=-8 pitch angle
   double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   
   double z0=parameters[6];                          
   double fz=exp(-fabs(z/z0));                                 
   if(z<0.) fz=-fz;                              // to avoid problem at z=0

   double  psi=1./(tan(p));//radial dependence of opening angle   

   double r=pow(x*x+y*y,0.5);
   double rbreak = parameters[8];                                  //EO20081119

    if(r<=rbreak) B0  =parameters[0]*(ro/rbreak);//parameters[0]=1.4
    if(r> rbreak) B0  =parameters[0]*(ro/r);

    Bregx=B0  *cos(theta-psi*log(r/r0))*sin(p-theta)*fz;                          
    Bregy=B0  *cos(theta-psi*log(r/r0))*cos(p-theta)*fz;
    Bregz=parameters[5];                                                              

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                                  

    // random component can be added for testing
    Bran=parameters[4];
    Branx=Bran*fz;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }


////////////////////////////////////////////////////////////////////////////////////////////
// halo dipole component from prouza 2003 formulae, same parameters

  if(name == "halo_dipole")
 {

        
    theta=atan2(y,x);             // NB theta defined differently from above
   double r=pow(x*x+y*y,0.5);
   double beta=90.*dtr-atan(z/r);

 if(fabs(z)<=0.3 && r<=0.1) //Btot=const =2mG from prouza
{double k=pow(10.,-3.);
Bregx=-(3./2.)*k*sin(2.*beta)*sin(theta);
Bregy=-(3./2.)*k*sin(2.*beta)*cos(theta);
Bregz=-k*(3.*cos(beta)*cos(beta)-1.);
}

 if(fabs(z)>0.3 && r<=0.1) //Btot=const =2mG from prouza 
{double k=2.*pow(10.,-7.);
Bregx=-(3./2.)*k*sin(2.*beta)*sin(theta)/pow(r,3.);
Bregy=-(3./2.)*k*sin(2.*beta)*cos(theta)/pow(r,3.);
Bregz=-k*(3.*cos(beta)*cos(beta)-1.)/pow(r,3.);
}

if( r>0.1 && r<=2.) // from prouza
{double k=2.*pow(10.,-7.);
Bregx=-(3./2.)*k*sin(2.*beta)*sin(theta)/pow(r,3.);
Bregy=-(3./2.)*k*sin(2.*beta)*cos(theta)/pow(r,3.);
Bregz=-k*(3.*cos(beta)*cos(beta)-1.)/pow(r,3.);
}

if( r>2. && r<=5.) // from prouza
{double k=5.*pow(10.,-7.);
Bregx=-(3./2.)*k*sin(2.*beta)*sin(theta);
Bregy=-(3./2.)*k*sin(2.*beta)*cos(theta);
Bregz=-k*(3.*cos(beta)*cos(beta)-1.);
}

if(  r>5.) // from prouza
{double k=pow(10.,-4.);
Bregx=-(3./2.)*k*sin(2.*beta)*sin(theta)/pow(r,3.);
Bregy=-(3./2.)*k*sin(2.*beta)*cos(theta)/pow(r,3.);
Bregz=-k*(3.*cos(beta)*cos(beta)-1.)/pow(r,3.);
}



    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                              

// random component set to zero
   Bran =0.0;
   Branx=0.0;  
   Brany=0.0;
   Branz=0.0;


}


/////////////////////////////////////////////////////// EO20090701
  //Sun model for B in the plane. ASS-RING model

 if(name == "sun")            //EO20090207
 {
   double B0;                                                                             

   double  ro=8.5;
   theta=atan2(y,x);
   double p=parameters[2]*dtr;                      //=-12 pitch angle
   //double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   double r0=10.;   

   double z0=parameters[6];         //=1 kpc                 
   double fz=exp(-fabs(z/z0));                                 
   if(z<0.) fz=-fz;                              // to avoid problem at z=0

   // double  psi=1./(tan(p));//radial dependence of opening angle   

   double r=pow(x*x+y*y,0.5);
   double rbreak = parameters[8];//rc=5 kpc                                  //EO20081119
   double D1;
   double D2;
   B0=parameters[0];//2 microGauss

   if(r<=rbreak && fabs(z)<=1) D1  = B0;
   if(r<=rbreak && fabs(z)>1) D1  = 0.;//artificial cut not present in Sun formulation 
   if(r> rbreak) D1  =B0*exp(-(r-ro)/r0)*fz;

   if(r>7.5) D2=1;
 if(r<=7.5 & r>6.) D2=-1;
 if(r>5.&& r<=6) D2=1;
 if(r<=5) D2=-1;


    Bregx=D1*D2*sin(p-theta);                          
    Bregy=D1*D2*cos(p-theta);  
    Bregz=parameters[5];//=0 in Sun                                                              

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                                  

    // random component can be added for testing
    Bran=parameters[4];//3 microgauss
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }

/////////////////////////////////////////////////////// AWS20101108
  //Sun etal. A&A 477, 573 (2008) model for B in the plane. ASS-RING model
 // alternative formulae with more explicit angular relations. Based on Sun etal 2008 equation (7).
 // Sun etal model has sun at x=-10 in a RH system. To get the galprop system with sun on +x.
 // The transformation x -> -x, y -> -y is used to goto the Sun etal system  and compute B_R and B_phi in the Sun etal system.
 // This is not actually needed for ASS_RING which has no azumithal modulation, but will be needed for ASS_ARM or BSS
 // These radial and azimuthal fields are the not affected by this rotation, so they are also valid in the galprop system.
 // Hence Bx and By can then be computed from these in the galprop system.

if(name == "Sun_ASS_RING")                            //AWS20101108
 {
   double B0      =parameters[0];      // 2 microGauss      in Sun etal. 
   double R0_ran  =parameters[1];      // R-scale       for random           AWS20110223
   double p       =parameters[2]*dtr;  // pitch angle  =-12 in Sun etal.  
   double z0_ran  =parameters[3];      // z-scaleheight for random           AWS20110222  
   double Bran0   =parameters[4];      // 3 microgauss      in Sun etal.     AWS20110222
          Bregz   =parameters[5];      // 0                 in Sun etal                                                                           
   double z0      =parameters[6];      // z-scaleheight for regular B. 1 kpc in Sun etal. 
   //              parameters[7]       reserved for halo field added for all models
   double Rc     = parameters[8];      // 5 kpc             in Sun etal  
   double Bc     = parameters[9];      // 2 microgauss      in Sun etal      AWS20101215

   double Rsun=8.5;   // solar position
   double R0  =10.;   // radial scale length

   double phi=atan2(-y,-x); // in the Sun etal system : not needed in ASS_RING but will be for other models



               
   double fz=exp(-fabs(z/z0));                                 
   

     

   double R=pow(x*x+y*y,0.5);
                             
   double D1;
   double D2;


   if(R<=Rc     && fabs(z)<=1) D1  = Bc;// replaces B0 to correspond to Sun etal formula AWS20101215
   if(R<=Rc     && fabs(z) >1) D1  = 0.;//artificial cut not present in Sun formulation 
   if(R> Rc    )               D1  = B0*exp(-(R-Rsun)/R0)*fz;

   if(R >7.5)           D2=+1;
   if(R<=7.5&& R >6.0)  D2=-1;
   if(R >5.0&& R<=6.0)  D2=+1;
   if(R<=5.0)           D2=-1;

   // formula (6) of Sun etal. A&A 477, 573 (2008)
   double B_R   =  D1 * D2 * sin(p); // radial    field 
   double B_phi = -D1 * D2 * cos(p); // azimuthal field in direction of increasing azimuth angle theta 

   // in galprop system use theta to project onto axes
   theta=atan2( y, x); // in the galprop  system

   Bregx = B_R*cos(theta) - B_phi*sin(theta);
   Bregy = B_R*sin(theta) + B_phi*cos(theta);

                                                           

   Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                                  

   // random component 

   double fR_ran=exp(-      (R-Rsun)/R0_ran );                                          //AWS20110223
   double fz_ran=exp(  -fabs(     z /z0_ran));                                          //AWS20110222

   Branx=Bran0 * fR_ran * fz_ran;  // put all field in x-direction (has no significance) AWS20110223
   Brany=0.0;
   Branz=0.0;
   Bran = sqrt(Branx*Branx+Brany*Brany+Branz*Branz);                      //AWS20110222 needed for 2D !

     if(debug==1)
     {
       cout<<"B_field_3D_model details: "<<name<<" x="<<x<<" y="<<y<<" z="<<z<<" B0="<<B0<<" R="<<R<<" theta="<<theta<<" B_R="<<B_R<<" B_phi="<<B_phi<<" Bregx="<<Bregx<<" Bregy="<<Bregy<<endl;
}

 }

/////////////////////////////////////////////////////// AWS20110311
//Sun etal. A&A 477, 573 (2008) model for B in the plane. ASS-RING model
// alternative formulae with more explicit angular relations. Based on Sun etal 2008 equation (7).
// Sun etal model has sun at x=-10 in a RH system. To get the galprop system with sun on +x.
// The transformation x -> -x, y -> -y is used to goto the Sun etal system  and compute B_R and B_phi in the Sun etal system.
// This is not actually needed for ASS_RING which has no azumithal modulation, but will be needed for ASS_ARM or BSS
// These radial and azimuthal fields are the not affected by this rotation, so they are also valid in the galprop system.
// Hence Bx and By can then be computed from these in the galprop system.
// This model is based on Sun_ASS_RING but has a more logical arrangement of parameters, and some fixed at Sun etal values.

if(name == "Sun_ASS_RING_2")                            //AWS20110311
 {
   // regular B

   double p       =parameters[0]*dtr;  // pitch angle  =-12 deg                      in Sun etal. 
   double B0      =parameters[1];      // B(Rsun)        for regular b. 2 microGauss in Sun etal. 
   double R0      =parameters[2];      // R-scale length for regular B. 10 kpc       in Sun etal.
   double z0      =parameters[3];      // z-scale height for regular B.  1 kpc       in Sun etal. 
 
         
   // random B 

   double Bran0   =parameters[4];      // B(Rsun)       for random B. 3 microgauss   in Sun etal.   
   double R0_ran  =parameters[5];      // R-scale       for random B.     infinite   in Sun etal.
   double z0_ran  =parameters[6];      // z-scaleheight for random B.     infinite   in Sun etal.
                                                                   

   //              parameters[7]       // reserved for halo field added for all models

   double Rc     = parameters[8];      // 5 kpc             in Sun etal  
   double Bc     = parameters[9];      // 2 microgauss      in Sun etal    

   // regular B fixed parameters  now use parameters AWS20110907
   /*
   double Bc      =          2.0e-6;   // 2 microgauss                                in Sun etal.      
   double Rc      =          5.0;      // 5 kpc                                       in Sun etal. 
   */

          Bregz   =parameters[10];                // 0                                in Sun etal. AWS20110907

   double Rsun=8.5;   // solar position
  
   double phi=atan2(-y,-x); // in the Sun etal system : not needed in ASS_RING but will be for other models

   double fz=exp(-fabs(z/z0));                                 
   
   double R=pow(x*x+y*y,0.5);
                             
   double D1;
   double D2;


   if(R<=Rc     && fabs(z)<=1) D1  = Bc;// replaces B0 to correspond to Sun etal formula AWS20101215
   if(R<=Rc     && fabs(z) >1) D1  = 0.;//artificial cut not present in Sun formulation 
   if(R> Rc    )               D1  = B0*exp(-(R-Rsun)/R0)*fz;

   if(R >7.5)           D2=+1;
   if(R<=7.5&& R >6.0)  D2=-1;
   if(R >5.0&& R<=6.0)  D2=+1;
   if(R<=5.0)           D2=-1;

   // formula (6) of Sun etal. A&A 477, 573 (2008)
   double B_R   =  D1 * D2 * sin(p); // radial    field 
   double B_phi = -D1 * D2 * cos(p); // azimuthal field in direction of increasing azimuth angle theta 

   // in galprop system use theta to project onto axes
   theta=atan2( y, x); // in the galprop  system

   Bregx = B_R*cos(theta) - B_phi*sin(theta);
   Bregy = B_R*sin(theta) + B_phi*cos(theta);

                                                           

   Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                                  

   // random component 

   double fR_ran=exp(-      (R-Rsun)/R0_ran );                                          //AWS20110223
   double fz_ran=exp(  -fabs(     z /z0_ran));                                          //AWS20110222

   Branx=Bran0 * fR_ran * fz_ran;  // put all field in x-direction (has no significance) AWS20110223
   Brany=0.0;
   Branz=0.0;
   Bran = sqrt(Branx*Branx+Brany*Brany+Branz*Branz);                      //AWS20110222 needed for 2D !

     if(debug==1)
     {
       cout<<"B_field_3D_model details: "<<name<<" x="<<x<<" y="<<y<<" z="<<z<<" B0="<<B0<<" R="<<R<<" theta="<<theta<<" B_R="<<B_R<<" B_phi="<<B_phi<<" Bregx="<<Bregx<<" Bregy="<<Bregy<<" Bregz="<<Bregz  <<endl;
}

 }



  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 // from Elena Orlando                                                                    EO20090716
  if(name == "random")
  {
 


    double Bo,rscale,zscale;
    double ro=8.5;
    double b_field;

    

    double r=sqrt(x*x+y*y);


   //just to test the random compoennt since the contrubution of the models for the regular is negligible EO20090716
   Bo   = parameters[0]; // Gauss
   rscale=parameters[3]; // kpc
   zscale=parameters[2]; // kpc

   double chi0=parameters[1]*dtr;// not defined in Han; parameter[1]=8 in WMAP Miville
    double chi=chi0*tanh(z/parameters[9]);      //  z dependence, not defined in Han;  //EO20081120

    b_field=Bo *exp(-(r-ro)/rscale)*cos(chi);// * exp(-fabs(z)/zscale);

 

   // the regular field is zero
   Breg =0.0;
   Bregx=0.0;
   Bregy=0.0;
   Bregz=0.0;

   Bran=b_field;
   Branx=Bran;  // put all field in x-direction (has no significance)
   Brany=0.0;
   Branz=0.0;

   theta=0.0;
   phi  =0.0;

  } //only random component

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 // from Elena Orlando                                                                    EO20090716
  if(name == "random1")
  {
 


    double Bo,rscale,zscale;
    double ro=8.5;
    double b_field;

    

    double r=sqrt(x*x+y*y);


   //just to test the random compoennt since the contrubution of the models for the regular is negligible EO20090716
   Bo   = parameters[0]; // Gauss
   rscale=parameters[3]; // kpc
   zscale=parameters[2]; // kpc

   double chi0=parameters[1]*dtr;// not defined in Han; parameter[1]=8 in WMAP Miville
    double chi=chi0*tanh(z/parameters[9]);      //  z dependence, not defined in Han;  //EO20081120
    double rc=parameters[4];

    if(r<=rc) b_field=Bo *exp(-(r-ro)/rscale)*cos(chi);// * exp(-fabs(z)/zscale);
    if(r>rc) b_field=Bo *exp(-(rc-ro)/rscale)*cos(chi);// * exp(-fabs(z)/zscale);
 

   // the regular field is zero
   Breg =0.0;
   Bregx=0.0;
   Bregy=0.0;
   Bregz=0.0;

   Bran=b_field;
   Branx=Bran;  // put all field in x-direction (has no significance)
   Brany=0.0;
   Branz=0.0;

   theta=0.0;
   phi  =0.0;

  } //only random component


/////////////////////////////////////////////////////// EO20090701
  //Sun model for B in the plane. ASS-RING model

 if(name == "jansson_disk")            //EO20090207
 {
   double B0;                                                                             

   double  ro=8.5;
   theta=atan2(y,x);
   double p=-5*dtr;                      //=-12 pitch angle
   //double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   double r0=5.1;   

   double z0=parameters[6];         //=1 kpc                 
   double fz=exp(-fabs(z/z0));                                 
   if(z<0.) fz=-fz;                              // to avoid problem at z=0

   // double  psi=1./(tan(p));//radial dependence of opening angle   

   double r=pow(x*x+y*y,0.5);
   double rbreak = 5.7;//rc=5 kpc                                  //EO20081119
   double D1;
   double D2;
   double Bc=0.16e-6;
   B0=1.1e-6;//2 microGauss

   if(r<=rbreak) D1  = Bc*fz;
   //  if(r<=rbreak) D1  = 0.;//artificial cut not present in Sun formulation 
   if(r> rbreak) D1  =B0*exp(-(r-ro)/r0)*fz;

   if(r>7.5) D2=1;
 if(r<=7.5 & r>6.) D2=-1;
 if(r>5.&& r<=6) D2=1;
 if(r<=5) D2=-1;


    Bregx=D1*D2*sin(p-theta);                          
    Bregy=D1*D2*cos(p-theta);  
    Bregz=0.0;//=0 in Sun                                                              

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                                  

    // random component can be added for testing
    Bran=0.0;//3 microgauss
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }


/////////////////////////////////////////////////////// EO20090701
  //Sun model for B in the plane. ASS-RING model
 //to be finished!!!!
 if(name == "jansson_halo")            //EO20090207
 {
   double B0;                                                                             

   double  ro=8.5;
   theta=atan2(y,x);
   double p;                      //=-12 pitch angle
   //double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   double r0=28.;   

   double z0=parameters[6];         //=1 kpc                 
   double fz=exp(-fabs(z/z0));                                 
   if(z<0.) fz=-fz;                              // to avoid problem at z=0

   // double  psi=1./(tan(p));//radial dependence of opening angle   

   double r=pow(x*x+y*y,0.5);
   double rbreak = 8.72;//rc=5 kpc                                  //EO20081119
   // double D1;
   // double D2;
   // double Bc=0.16;
   B0=2.3e-6;//2 microGauss

   if(r<=rbreak) p=-30.*dtr;//     D1  = Bc*fz;
   //  if(r<=rbreak) D1  = 0.;//artificial cut not present in Sun formulation 
   if(r> rbreak) p=-2.*dtr;  // D1  =B0*exp(-r/r0)*fz;

   //   if(r>7.5) D2=1;
   // if(r<=7.5 & r>6.) D2=-1;
   // if(r>5.&& r<=6) D2=1;
   // if(r<=5) D2=-1;


    Bregx=B0*exp(-r/r0)*fz*sin(p-theta);                          
    Bregy=B0*exp(-r/r0)*fz*cos(p-theta);  
    Bregz=0.0;//=0 in Sun                                                              

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                                  

    // random component can be added for testing
    Bran=0.0;//3 microgauss
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }


//----------------------------------------------------------------------
if(name == "Pshirkov_ASS")                            //EO20110925
 
  // Pschirkov, Tinyakov, Kronberg & Newton-McGee,  2011 ApJ 738,192
  // axisymmetric spiral, no reversals in disk
 {
   // regular B

   double p       =parameters[0]*dtr;  // pitch angle  = -5 deg                           in Ps etal. 
   double B0      =parameters[1];      // B(Rsun)        for regular b. 2 microGauss      in Ps etal. 
   double d       =parameters[2];      // distance to the first field reversal -0.6 kpc   in Ps etal.
   double z0      =parameters[3];      // z-scale height for regular B.  1 kpc            in Ps etal. 
 
   double Rc     = parameters[8];      // B=B(Rc) R<Rc                   5 kpc            in Ps etal.  
   //            = parameters[9];      // not used                                        
         
   // random B 

   double Bran0   =parameters[4];      // B(Rsun)       for random B. 3 microgauss   in Sun etal.   
   double R0_ran  =parameters[5];      // R-scale       for random B.     infinite   in Sun etal.
   double z0_ran  =parameters[6];      // z-scaleheight for random B.     infinite   in Sun etal.
                                                                   

   //              parameters[7]       // reserved for halo field added for all models


          Bregz   =parameters[10];              

   double Rsun=8.5;   // solar position
  

   double fz=exp(-fabs(z/z0));                                 
   
   double R=pow(x*x+y*y,0.5);
        
                    
   double b  = 1./tan ( p );    
   double phi= b*log(1. + d/Rsun)-pi/2.;
   

   double    B_r; 
   if(R<=Rc) B_r  = B0*Rsun/(Rc*cos(phi)) ;
   if(R> Rc) B_r  = B0*Rsun/(R *cos(phi)) ;           

   // in galprop system use theta to project onto axes
   theta=atan2( y, x); // azimuth angle in the galprop  system 


   // equation (3) of Pshirkov et al.
   double theta_ps = -theta; // since Ps uses opposite convention for theta
   R += 1.0e-6; // avoid log(0) !
   double B=B_r* fabs( cos(theta_ps-b*log(R/Rsun)+phi)) *fz; // only difference from Pshirkov BSS 


   double B_R     =  B * sin(p); // radial    field 


   double B_theta = -B * cos(p); // azimuthal field in direction of increasing azimuth angle theta

   //         Ps has B * cos(p) but this seems because they define azimuth clockwise, while we have anticlockwise.
   // see Tinyakov 2002 APh 18,165: "local field points to l=90+p" so p=-5 deg gives l=85 and hence clockwise from above.
   // so to get local B clockwise in our system, need minus (like Sun etal).
   // Ps base their system on Han and Qiao 1994 A&A 288,759 which has a diagram with azimuth clockwise, hence confirmed.



   Bregx = B_R*cos(theta) - B_theta*sin(theta);
   Bregy = B_R*sin(theta) + B_theta*cos(theta);

                                                           

   Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                                  

   // random component 

   double fR_ran=exp(-      (R-Rsun)/R0_ran );                                         
   double fz_ran=exp(  -fabs(     z /z0_ran));                                          

   Branx=Bran0 * fR_ran * fz_ran;  // put all field in x-direction (has no significance)
   Brany=0.0;
   Branz=0.0;
   Bran = sqrt(Branx*Branx+Brany*Brany+Branz*Branz);                      //            needed for 2D !


   if(debug==2)// not yet invoked via  galdef verbose
   {
     cout<<"B_field_3D_model  name= "<<name
     <<" (x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<" theta= "<<theta
     <<"  B_r ="<< B_r <<"  B ="<< B <<" B_R  ="<<  B_R  <<" B_theta="<< B_theta<<endl
     <<"  theta_ps ="<< theta_ps <<"  b ="<< b <<" phi  ="<<  phi 
     <<" theta_ps-b*log(R/Rsun)+phi ="<<  theta_ps-b*log(R/Rsun)+phi
     <<" cos(theta_ps-b*log(R/Rsun)+phi) ="<<  cos(theta_ps-b*log(R/Rsun)+phi)<<endl
     <<"  Breg="<< Breg <<" Bregx="<<  Bregx<<" Bregy="<< Bregy<<" Bregz="<< Bregz
     <<" Bran="<< Bran <<" Branx="<<  Branx<<" Brany="<< Brany<<" Branz="<< Branz<<endl<<endl;
   }
 }


   //------------------------------------------------------------------------------------------------------------
if(name == "Pshirkov_BSS")                            //EO20110925

  // Pschirkov, Tinyakov, Kronberg & Newton-McGee,  2011 ApJ 738,192
  // bisymmetric spiral, with reversals in disk
 {
   // regular B

   double p       =parameters[0]*dtr;  // pitch angle  = -6 deg                           in Ps etal. 
   double B0      =parameters[1];      // B(Rsun)        for regular b. 2 microGauss      in Ps etal. 
   double d       =parameters[2];      // distance to the first field reversal -0.6 kpc   in Ps etal.
   double z0      =parameters[3];      // z-scale height for regular B.  1 kpc            in Ps etal.   
   double Rc     = parameters[8];      // B=B(Rc) R<Rc                   5 kpc            in Ps etal.  
   //            = parameters[9];      // not used                                      
         
   // random B 

   double Bran0   =parameters[4];      // B(Rsun)       for random B. 3 microgauss   in Sun etal.   
   double R0_ran  =parameters[5];      // R-scale       for random B.     infinite   in Sun etal.
   double z0_ran  =parameters[6];      // z-scaleheight for random B.     infinite   in Sun etal.
                                                                   

   //              parameters[7]       // reserved for halo field added for all models


          Bregz   =parameters[10];                

   double Rsun=8.5;   // solar position
  
   double fz=exp(-fabs(z/z0));                                 
   
   double R=pow(x*x+y*y,0.5);
        
   double b  = 1./tan ( p );   
   double phi= b*log(1. + d/Rsun)-pi/2.;
   

   double    B_r;    
   if(R<=Rc) B_r  = B0*Rsun/(Rc*cos(phi)) ;
   if(R> Rc) B_r  = B0*Rsun/(R *cos(phi)) ;           



   // in galprop system use theta to project onto axes
   theta=atan2( y, x); // azimuth angle in the galprop  system 

  

   // equation (4) of Pshirkov etal.
   double theta_ps = -theta; // since Ps uses opposite convention for theta
   R += 1.0e-6; // avoid log(0) !
   double B=B_r*  cos(theta_ps-b*log(R/Rsun)+phi) * fz; // only difference from Pshirkov ASS 


     // formula (3) of Ps et al.
   double B_R     =  B * sin(p); // radial    field 
   double B_theta = -B * cos(p); // azimuthal field in direction of increasing azimuth angle theta 

   //         Ps has B * cos(p) but this seems because they define azimuth clockwise, while we have anticlockwise.
   // see Tinyakov 2002 APh 18,165: "local field points to l=90+p" so p=-5 deg gives l=85 and hence clockwise from above.
   // so to get local B clockwise in our system, need minus (like Sun etal).
   // Ps base their system on Han and Qiao 1994 A&A 288,759 which has a diagram with azimuth clockwise, hence confirmed.



   Bregx = B_R*cos(theta) - B_theta*sin(theta);
   Bregy = B_R*sin(theta) + B_theta*cos(theta);

                                                           

   Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                                  

   // random component 

   double fR_ran=exp(-      (R-Rsun)/R0_ran );                                          
   double fz_ran=exp(  -fabs(     z /z0_ran));                                         

   Branx=Bran0 * fR_ran * fz_ran;  // put all field in x-direction (has no significance) 
   Brany=0.0;
   Branz=0.0;
   Bran = sqrt(Branx*Branx+Brany*Brany+Branz*Branz);                      //            needed for 2D !


   if(debug==2) // not yet invoked via  galdef verbose
   {
     cout<<"B_field_3D_model  name= "<<name
     <<" (x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<" theta= "<<theta
     <<"  B_r ="<< B_r <<"  B ="<< B <<" B_R  ="<<  B_R  <<" B_theta="<< B_theta<<endl
     <<"  theta_ps ="<< theta_ps <<"  b ="<< b <<" phi  ="<<  phi 
     <<" theta_ps-b*log(R/Rsun)+phi ="<<  theta_ps-b*log(R/Rsun)+phi
     <<" cos(theta_ps-b*log(R/Rsun)+phi) ="<<  cos(theta_ps-b*log(R/Rsun)+phi)<<endl
     <<"  Breg="<< Breg <<" Bregx="<<  Bregx<<" Bregy="<< Bregy<<" Bregz="<< Bregz
     <<" Bran="<< Bran <<" Branx="<<  Branx<<" Brany="<< Brany<<" Branz="<< Branz<<endl<<endl;
   }

 
}

   //------------------------------------------------------------------------------------------


double Bhalox,Bhaloy,Bhaloz; //AWS20110926

if(name != "Pshirkov_ASS" && name != "Pshirkov_BSS") //AWS20110926
{
 /////////////////////////////////////////////////////////////////////////////////////////////
 // from Elena Orlando                                                                    EO20081114
 //  Sun et al. 2008, A&A 477, 573; Prouza & Smida 2003 A&A 410, 1
 //  updated to Sun & Reich, arxiv:1010.4394 by AWS20101110
 //  halo field, added to regular field

 // always added but controlled by parameters[7] which can be set to zero if halo not required
 // so this code should come at end

 double B0t=parameters[7]; //2 microgauss in Sun & Reich, arxiv:1010.4394

 double r=pow(x*x+y*y,0.5);

 double fi=atan2(r,z);
 
 double z0 =1.5;                                      //AWS20101110
 double R0 =4.0;                                      //AWS20101110
 double z1 =0.;
 if (fabs(z) <1.5) z1=0.2;
 if (fabs(z)>=1.5) z1=4.0;//  updated to Sun & Reich, arxiv:1010.4394 by AWS20101110

// if (fabs(z)>=1.5) z1=0.4; Sun et al. 2008, A&A 477, 573, now superseded

 double B1=1.0/(1.0 + pow((fabs(z)-z0)/z1, 2));             //AWS20101110

 double B2=exp(-(r-R0)/R0);                                 //AWS20101110
 double Bt=B0t*B1*B2*r/R0;                                  //AWS20101110

 theta=atan2(y,x); // azimuth angle in galprop system



// B is defined as in increasing theta direction (as for Sun_ASS_RING)
// B0t should then be +ve to get anticlockwise field for z>0 as required

        Bhalox = -Bt*sin(theta);                                 //AWS20110926
        Bhaloy =  Bt*cos(theta);                                 //AWS20110926

// B reversal below plane
 if (z<0.){Bhalox=-Bhalox; Bhaloy=-Bhaloy;}

        Bhaloz=0.;



 } //name not Pshirkov_ASS, BSS

 /////////////////////////////////////////////////////////////////////////////////////////////



if(name == "Pshirkov_ASS" || name == "Pshirkov_BSS") //AWS20110926
{
 // from Elena Orlando HALO FIELD                                                                   EO20110926
// updated Ps et al. ApJ 738,192
 //  halo field, added to regular field

 // always added but controlled by parameters[7] which can be set to zero if halo not required
 // so this code should come at end

 double BH0=parameters[7]; //4 microgauss in Ps et al. ApJ 738,192

 double r=pow(x*x+y*y,0.5);

 double fi=atan2(r,z);
 
 double zH0 =1.3;                                   // Ps et al. ApJ 738,192  Table 3
 double RH0 =8.0;                                   // Ps et al. ApJ 738,192  Table 3
 double zH1 =0.;

 if (fabs(z)< zH0 ) zH1=0.25;                       // Ps et al. ApJ 738,192  Table 3
 if (fabs(z)>=zH0 ) zH1=0.4;                        // Ps et al. ApJ 738,192  Table 3

 // only exception: Ps etal Table 3.  halo South, ASS  BH0=2 when North=4: apply same factor to user parameter
 if(name == "Pshirkov_ASS" && z<0.0) BH0/=2.0;



 double B1=1.0/(1.0 + pow((fabs(z)-zH0)/zH1, 2));            
 double B2=exp(-(r-RH0)/RH0);                                
 double Bt=BH0 *B1*B2*r/RH0;                        // Ps et al. ApJ 738,192 equation (8)                                                        

 theta=atan2(y,x); // azimuth angle in galprop system



// B is defined as in increasing theta direction 
// BH0 should then be +ve to get anticlockwise field for z>0 as required

        Bhalox = -Bt*sin(theta);                                 //AWS20110926
        Bhaloy =  Bt*cos(theta);                                 //AWS20110926

// B reversal below plane

// if (z<0.){Bhalox=-Bhalox/2.; Bhaloy=-Bhaloy/2.;} // the field in the south is 1/2!!!!
   if (z<0.){Bhalox=-Bhalox;    Bhaloy=-Bhaloy   ;} //AWS20120405 1/2 factor already applied above, applies to ASS only

   Bhaloz=0.;




 } //name is  Pshirkov_ASS, BSS

// this is always required:


 Bregx += Bhalox;
 Bregy += Bhaloy;
 Bregz += Bhaloz;

 Breg=sqrt(Bregx*Bregx + Bregy*Bregy +  Bregz*Bregz); 

//////////////////////////////////////////////////////////////// 

 


///////////////////////////////////////////////////////AWS20140516  (from M.Regis and M.Taoso)
// adapted from Fornengo etal. http://arxiv.org/abs/1402.2218 routine available at http://www.astroparticle.to.infn.it/darkmatter/Research/Radio.1402.2218/

  //Farrar and Jansson 2012 model for regular and random B (see 1204.3662 and 1210.7820)
 if(name == "JF12_Fornengo")      //AWS20140516 the only change from the original: name is string instead of char      
 {
  double bregdisk[3],bregtor[3],bregXhalo[3],brandisk,branhalo;
  double i=11.5*pi/180.;
  double bring=0.1,hdisk=0.4,wdisk=0.27;	
  double b1=0.1,b2=3.0,b3=-0.9,b4=-0.8,b5=-2.0,b6=-4.2,b7=0.0,b8=2.7;
  double Rx1=-5.1,Rx2=-6.3,Rx3=-7.1,Rx4=-8.3,Rx5=-9.8,Rx6=-11.4,Rx7=-12.7,Rx8=-15.5;
  double r,rsph,phi2,B0,bt,rt,phil;
  double Lf=1.0/(1.0+exp(-2*(fabs(z)-hdisk)/wdisk));
  if(x!=0.)	{phi2=atan(fabs(y/x));}
  r=sqrt(x*x+y*y);
  rsph=sqrt(x*x+y*y+z*z);
  if(x==0. && y>=0.0)	{phi=pi/2.;}
  if(x==0. && y<0.0)	{phi=3.*pi/2.;}
  if(x>0. && y>=0.0)	{phi=phi2;}
  if(x<0. && y>=0.0)	{phi=pi-phi2;}
  if(x<0. && y<0.0)	{phi=pi+phi2;}
  if(x>0. && y<0.0)	{phi=2*pi-phi2;}

// Regular B: disk component:
  if(r>=20.0)	{B0=0;}
  if(r<3.0)	{B0=0;}
  if(r>=3.0 && r<=5.0)	{B0=bring;}
  if(r>5. && r<20.){

    if(phi<=pi)	{phil=pi+phi;} else {phil=phi-pi;}
    rt=r*exp(-phil/(tan(pi/2-i)));
    if(rt<fabs(Rx1))	{rt=r*exp(-(phil-2*pi)/tan(pi/2.-i));}
    if(rt>fabs(Rx8))	{rt=r*exp(-(phil+2*pi)/tan(pi/2.-i));}

    if(rt<fabs(Rx1))	bt=b1;
    if((rt>fabs(Rx1))&&(rt<fabs(Rx2)))	bt=b2;
    if((rt>fabs(Rx2))&&(rt<fabs(Rx3)))	bt=b3;
    if((rt>fabs(Rx3))&&(rt<fabs(Rx4)))	bt=b4;
    if((rt>fabs(Rx4))&&(rt<fabs(Rx5)))	bt=b5;
    if((rt>fabs(Rx5))&&(rt<fabs(Rx6)))	bt=b6;
    if((rt>fabs(Rx6))&&(rt<fabs(Rx7)))	bt=b7;
    if((rt>fabs(Rx7))&&(rt<fabs(Rx8)))	bt=b8;

  B0=bt*5.0/r*(1-Lf);
  }

  double br=0.0;
  double bphi=1.0;
  if(r>5.){
    br=sin(i);
    bphi=cos(i);}	

   bregdisk[0]=(br*cos(phi)-cos(pi/2.-phi)*bphi)*B0;
   bregdisk[1]=(br*sin(phi)+sin(pi/2.-phi)*bphi)*B0;
   bregdisk[2]=0.0;


// Regular B: toroidal component:
  double bn=1.4,bs=-1.1,rn=9.22,rs=16.7,wh=0.2,z0=5.3;
  double Lfn=1.0/(1.+exp(-2*(fabs(r)-rn)/wh));
  double Lfs=1.0/(1.+exp(-2*(fabs(r)-rs)/wh));
  B0=(z<0) ? exp(-fabs(z/z0))*Lf*bs*(1-Lfs) : exp(-fabs(z/z0))*Lf*bn*(1-Lfn); 
  if(rsph<1.0 || rsph>20.)	{B0=0;}

  br=0.0;
  bphi=1.0;
  bregtor[0]=(br*cos(phi)-cos(pi/2.-phi)*bphi)*B0;
  bregtor[1]=(br*sin(phi)+sin(pi/2.-phi)*bphi)*B0;
  bregtor[2]=0;

// Regular B: X-halo component:
  double thetaxu,theta0x=49.*pi/180.,rcx=4.8;	
  double bx=4.6,rx=2.9;
  double rp=r-fabs(z)/tan(theta0x);
  double BXp=bx*exp(-rp/rx);
  if(r==0)	{B0=0;thetaxu=0.0;}
  else	{B0=BXp*rp/r;thetaxu=theta0x;}

  if(rp<rcx){
    rp=r*rcx/(rcx+fabs(z)/tan(theta0x));
    if(r==rp )	{thetaxu=pi/2.;}
    else	{thetaxu=atan(fabs(z)/(r-rp));}
    B0=BXp*pow(rcx/(rcx+fabs(z)/tan(theta0x)) ,2.0);
  }
  if(rsph<1.0 || rsph>20.) B0=0;

  br=1.0; bphi=0.0;
  bregXhalo[0]=(br*cos(phi)-cos(pi/2.-phi)*bphi)*B0*cos(thetaxu);
  if(z<0) bregXhalo[0]*=-1;
  bregXhalo[1]=(br*sin(phi)+sin(pi/2.-phi)*bphi)*B0*cos(thetaxu);
  if(z<0) bregXhalo[1]*=-1;
  bregXhalo[2]=B0*sin(thetaxu);

// Regular B: Total
  Bregx=bregdisk[0]+bregtor[0]+bregXhalo[0];                          
  Bregy=bregdisk[1]+bregtor[1]+bregXhalo[1]; 
  Bregz=bregdisk[2]+bregtor[2]+bregXhalo[2];                                                              
  Breg =sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                                  

// random striated component added as a multiplicative factor acting on Breg:
  double fstri=sqrt(2.36); // the value is taken from the most recent paper 1210.7820
  Bregx*=fstri*1.e-6;
  Bregy*=fstri*1.e-6;
  Bregz*=fstri*1.e-6;
  Breg *=fstri*1.e-6; // Gauss


// Random B: disk component:
  double bint=7.63,z0disk=0.61;
  if(r<5.0)	{bt=bint;}
  if(r>5.){
    b1=10.81;b2=6.96;b3=9.59;b4=6.96;b5=1.96;b6=16.34;b7=37.29;b8=10.35;
    if(phi<=pi)	{phil=pi+phi;} else {phil=phi-pi;}
    rt=r*exp(-phil/(tan(pi/2-i)));
    if(rt<fabs(Rx1)) rt=r*exp(-(phil-2*pi)/tan(pi/2.-i));
    if(rt>fabs(Rx8)) rt=r*exp(-(phil+2*pi)/tan(pi/2.-i));
    if(rt<fabs(Rx1))	bt=b1;
    if((rt>fabs(Rx1))&&(rt<fabs(Rx2)))	bt=b2;
    if((rt>fabs(Rx2))&&(rt<fabs(Rx3)))	bt=b3;
    if((rt>fabs(Rx3))&&(rt<fabs(Rx4)))	bt=b4;
    if((rt>fabs(Rx4))&&(rt<fabs(Rx5)))	bt=b5;
    if((rt>fabs(Rx5))&&(rt<fabs(Rx6)))	bt=b6;
    if((rt>fabs(Rx6))&&(rt<fabs(Rx7)))	bt=b7;
    if((rt>fabs(Rx7))&&(rt<fabs(Rx8)))	bt=b8;
    bt=bt*5.0/r;
  }
  brandisk=bt*exp(-z*z/(2*z0disk*z0disk));	


// Random B: halo component:
  double r0=10.97;
  z0=2.84;
  B0=4.68;
  branhalo=B0*exp(-r/r0)*exp(-z*z/(2*z0*z0));

// Random B: Total (assumed to be isotropic)
  Bran=(brandisk+branhalo)*1.e-6; // Gauss
//  Bran=parameters[0]*exp(-(r-8.5)/parameters[1])*exp(-fabs(z)/parameters[2]) ; // uncomment this if you want to set the random component with a double-exponential law !
  Branx=Bran; // dummy! 
  Brany=0;
  Branz=0;


//  cout<<x<<" "<<y<<" "<<z<<" Btot Farrar "<<sqrt(Breg*Breg + Bran*Bran)<<endl;
// debug=2;

 if(debug==1)
 {
  cout<<"B_field_3D_mode name= "<<name
     <<" (x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<" theta= "<<theta                <<endl
     <<" (x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<" bregdisk  ="<<  bregdisk [0] <<"  "<< bregdisk [1]<<" "<< bregdisk [2]  <<" microG"     <<endl
     <<" (x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<" bregtor   ="<<  bregtor  [0] <<"  "<< bregtor  [1]<<" "<< bregtor  [2]  <<" microG"     <<endl
     <<" (x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<" bregXhalo ="<<  bregXhalo[0] <<"  "<< bregXhalo[1]<<" "<< bregXhalo[2]  <<" microG"     <<endl
     <<" (x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<" Breg="<< Breg <<" Bregx="<<  Bregx<<" Bregy="<< Bregy<<" Bregz="<< Bregz<<"      G"     <<endl
     <<" (x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<" brandisk ="<<  brandisk <<" branhalo ="<< branhalo                      <<" microG"     <<endl
     <<" (x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<" Bran ="<<  Bran <<" Branx ="<< Branx<<" Brany ="<< Brany<<" Branz ="<< Branz <<" G"     <<endl;
 }

} // farrar


//-----------------------------------------

 if(name == "JF12_original")      //AWS20140603
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

} //JF12

//============================================================================





if(debug==1)
  {
 cout<<"B_field_3D_model including halo: name= "<<name
     <<" (x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<" theta= "<<theta
     <<"  Breg="<< Breg <<" Bregx="<<  Bregx<<" Bregy="<< Bregy<<" Bregz="<< Bregz
                         <<" Bran="<< Bran <<" Branx="<<  Branx<<" Brany="<< Brany<<" Branz="<< Branz<<endl;
 
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
double Bperp(double x,double y,double z, double Bx,double By,double Bz, double x0,double y0,double z0)
{
  //returns the field (at x,y,z)  component perpendicular to the line-of sight to x0,y0,z0

  double Bperp_,Bperp2;

  // component parallel to line of sight by projection onto line of sight
  double Bpar=        ((x-x0)*Bx +     (y-y0)*By +    (z-z0)*Bz    ) 
                /sqrt ((x-x0)*(x-x0) + (y-y0)*(y-y0)+ (z-z0)*(z-z0)) ;

  // component perpendicular by Pythagoras: sqrt(B_total^2 - Bpar^2)
  double B_total2=(Bx*Bx + By*By + Bz*Bz);

  Bperp2 =      Bx*Bx + By*By + Bz*Bz - Bpar*Bpar ;
  Bperp_= 0.0;
  if(Bperp2 > 0.0)          // to avoid rounding error problem
  Bperp_ = sqrt(Bperp2);

  //  cout<<"Bperp routine: B_total^2="<<B_total2<< " Bpar="<<Bpar <<" Bperp2="<<Bperp2<<" Bperp_="<<Bperp_<<endl;

  

  return Bperp_;
}


//////////////////////////////////////////////////////////////////////////////////////////////

double B_field_3D_model_tot
(const std::string &name, const std::vector<double> &parameters,
 double x, double y, double z,
 int debug=0 )
{
  // return just total B, used e.g. for energy losses

  int options=0;
  double Breg,   Bregx,  Bregy,  Bregz;
  double Bran,   Branx,  Brany,  Branz;
  double B_tot;


  //  debug=1; // test AWS20080331

 B_field_3D_model
( name, parameters,
  x,  y,  z,
  options,
  Breg,  Bregx, Bregy, Bregz,
  Bran,  Branx, Brany, Branz,
  debug );

 B_tot=sqrt(Breg*Breg + Bran*Bran);

 return B_tot;

}

//////////////////////////////////////////////////////////////////////////////////////////////

double B_field_3D_model_tot
(const std::string &name, const std::vector<double> &parameters,
 double R, double z,
 int debug=0 )
{
  // return just total B, used e.g. for energy losses

  int options=0;
  double x,y;
  double Breg,   Bregx,  Bregy,  Bregz;
  double Bran,   Branx,  Brany,  Branz;
  double B_tot;

  // compute field on GC-Sun axis since x,y undefined (better would be to take azimuthal average)
  x = R;
  y = 0.0;

 B_field_3D_model
( name, parameters,
  x,  y,  z,
  options,
  Breg,  Bregx, Bregy, Bregz,
  Bran,  Branx, Brany, Branz,
  debug );

 B_tot=sqrt(Breg*Breg + Bran*Bran);

 if (debug==1) cout<<" B_field_3D_model_tot: R="<<R<<" z="<<z<<" Btot="<<B_tot<<endl; //AWS20101108

 return B_tot;

}

///////////////////////////////////////////////////////////////
//             test program
///////////////////////////////////////////////////////////////
/*

#include "synchrotron_emissivity.h"
int main()
{

  char name[20];
  double parameters[20];
  double x,y,z;
  int options;
  double  Breg,          Bregx,         Bregy,        Bregz;
  double  Bran,          Branx,         Brany,        Branz;
  int debug;

  strcpy(name,"test");
  strcpy(name,"circular");

  x=8;
  y=1;
  z=1;

  options=0;
  debug=1;

  for (x= -10.; x<+10.; x+=2) // avoid x0,y0,z0
  for (y= -10.; y<+10.; y+=2)
  for (z=  -1; z< +1.; z+=1.) 
  {

 B_field_3D_model
( name, parameters,
  x,  y,  z,
  options,
 Breg,  Bregx, Bregy, Bregz,
  Bran,  Branx, Brany, Branz,
  debug );

 cout<<"(x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  ";
 cout<<" B model= "<<name<<"  Breg="<<  Breg <<" Bregx="<<  Bregx<<" Bregy="<< Bregy<<" Bregz="<< Bregz
                         <<" Bran="<<  Bran <<" Branx="<<  Branx<<" Brany="<< Brany<<" Branz="<< Branz<<endl;

 double x0,y0,z0;
 x0=8.5; // solar position
 y0=0;
 z0=0;

 cout<<"Bperp for regular field ="<<Bperp( x, y, z,  Bregx, Bregy, Bregz,  x0,y0,z0)<<endl;
 
 cout<<"Bperp for random  field ="<<Bperp( x, y, z,  Branx, Brany, Branz,  x0,y0,z0)<<endl;

 // test synchrotron routine

 double gamma,nu;
 double synch_emissivity_total;
 double synch_emissivity_reg;
 double synch_emissivity_par;
 double synch_emissivity_perp;
 double synch_emissivity_random;

 int debug=1;

 gamma=1e4/.511; // 10 GeV electrons = 1e4 MeV
 nu=1.e9; // 1000 MHz


 double Brand; // NB names inconsistent
 Brand=Bran;
 double Bperp_reg=Bperp( x, y, z,  Bregx, Bregy, Bregz,  x0,y0,z0);
 double Bperp_ran=Bperp( x, y, z,  Branx, Brany, Branz,  x0,y0,z0);


 Bperp_reg += 1.e-9; // smaller field gives error in synchrotron.c
 Bperp_ran += 1.e-9; // smaller field gives error in synchrotron.c

 cout<<endl;
 cout<<"testing synchrotron routine"<<endl;
 cout<<endl;

 cout<<"regular field:"<<endl;
 cout<<"gamma= "<<gamma<<" nu="<<nu<<" Bperp_reg="<<Bperp_reg   << " Brand="<<Brand   <<endl;

 synch_emissivity_total
     = synchrotron_emissivity( gamma, nu, Bperp_reg, Brand,
                         synch_emissivity_reg, synch_emissivity_par, synch_emissivity_perp, synch_emissivity_random, debug );


 cout<<"regular field:"<<endl;
 cout<<"gamma= "<<gamma<<" nu="<<nu<<" Bperp_reg="<<Bperp_reg   << " Brand="<<Brand   <<endl

  <<" synch emissivity random   =  "  << synch_emissivity_random <<endl
  <<" synch emissivity reg      =  "  << synch_emissivity_reg    <<endl
  <<" synch emissivity parallel =  "  << synch_emissivity_par    <<endl
  <<" synch emissivity perp     =  "  << synch_emissivity_perp   <<endl 
  <<" synch emissivity total    =  "  << synch_emissivity_total  <<endl       
  <<endl;


 cout<<"random  field:"<<endl;
 cout<<"gamma= "<<gamma<<" nu="<<nu<<" Bperp_ran="<<Bperp_ran   << " Brand="<<Brand   <<endl;

 synch_emissivity_total
     = synchrotron_emissivity( gamma, nu, Bperp_ran, Brand,
                         synch_emissivity_reg, synch_emissivity_par, synch_emissivity_perp, synch_emissivity_random, debug );


 cout<<endl;
 cout<<"random  field:"<<endl;
 cout<<"gamma= "<<gamma<<" nu="<<nu<<" Bperp_ran="<<Bperp_ran   << " Brand="<<Brand   <<endl

  <<" synch emissivity random   = "  << synch_emissivity_random <<endl
  <<" synch emissivity reg      = "  << synch_emissivity_reg    <<endl
  <<" synch emissivity parallel = "  << synch_emissivity_par    <<endl
  <<" synch emissivity perp     = "  << synch_emissivity_perp   <<endl 
  <<" synch emissivity total    = "  << synch_emissivity_total  <<endl       
  <<endl;


  }// for x y z

 return 0;

}
*/
