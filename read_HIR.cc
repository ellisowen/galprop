
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * read_HIR.cc *                                galprop package * 5/22/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"
#include"fitsio.h"
#include"Skymap.h"
#include <vector>
#include <algorithm>
#include"SkymapFitsio.h"

int Galprop::read_HIR()
{
   cout<<" >>>> read_HIR    "<<endl;
   int status=0;

   fitsfile *fptr;
   char HIR_filename[200];

	 /*
	 if (galdef.skymap_format == 3) {
		 //Every ring stored in a separate file of order 9 (Nside = 512)
		 //The scheme is nested so there is an easy way to convert to maps of
		 //different sizes
		 const int order = 9;
		 const int n_Ring = 9; //Number of rings
		 double nulval = 0;
		 int anynul = 0;
		 int hdutype;
		 int nrow=1024; //Number of elements in each row in the table
		 //Resize the skymap and use NESTED scheme
	   galaxy.hpHIR.Resize(galdef.healpix_order, n_Ring, NEST);
		 //Number of pixels in read in map that is equivalent to one pixel in the
		 //used map
		 int nside = 1<<order;
		 long npix = 12*nside*nside;
		 const int pixelRatio = npix/galaxy.hpHIR.Npix();
		 const int invPixelRatio = galaxy.hpHIR.Npix()/npix;
		 //The fits directory
		 char HIR_file[200];
		 for (int i_Ring=0; i_Ring < n_Ring; ++i_Ring){
			 sprintf(HIR_file, "rbands_hi7_%d_512.fits", i_Ring);
			 strcpy(HIR_filename,configure.fits_directory);
			 strcat(HIR_filename,HIR_file);
			 if( fits_open_file(&fptr,HIR_filename,READONLY,&status) ) {
				 cout<<"read HIR open status= "<<status<<endl;
				 cout<<"HIR filename= "<<HIR_filename<<endl;
				 return(status);
			 }
			 //Move to the correct HDU number 2
			 if( fits_movabs_hdu(fptr,2,&hdutype,&status) ) {
				 cout<<"move to HIR binary table status= "<<status<<endl;
				 return status;
			 }
			 //How we read the map depends on which map has the finer grid.
			 //pixelRatio will be zero if hpHIR has finer grid
			 if (pixelRatio > 1){
				 double *tmparray = new double[pixelRatio];
				 //Loop through the hpHIR map and add it one pixel at a time
				 for( int p=0; p<galaxy.hpHIR.Npix(); ++p ){
					 long frow = p*pixelRatio/nrow + 1;
					 long felement = p*pixelRatio%nrow + 1;
					 fits_read_col(fptr,TDOUBLE,1,frow,felement,pixelRatio,&nulval,tmparray,&anynul,&status);
					 if (status){
						 cout<<"read HIR read table failed, status= "<<status<<endl;
						 cout<<"p= "<<p<<", pixelRatio= "<<pixelRatio<<", felement= "<<felement<<endl;
						 return(status);
					 }
					 //Calculate the average and insert into skymap
					 double sum = 0;
					 for( int i=0; i<pixelRatio; ++i){
						 sum += tmparray[i];
					 }
					 galaxy.hpHIR[p][i_Ring] = sum/double(pixelRatio);
				 }
				 delete[] tmparray;
			 }else{
				 double *tmparray = new double[npix];
				 //Read the whole skymap
				 fits_read_col(fptr,TDOUBLE,1,1L,1L,npix,&nulval,tmparray,&anynul,&status);
				 if (status){
					 cout<<"read HIR read table failed, status= "<<status<<endl;
					 return(status);
				 }
				 //Loop through the array and insert into 
				 for(int i=0; i<npix; ++i){
					 for(int p=i*invPixelRatio; p<(i+1)*invPixelRatio; ++p){
						 galaxy.hpHIR[p][i_Ring] = tmparray[i];
					 }
				 }
				 delete[] tmparray;
			 }
			 fits_close_file(fptr, &status);
			 if(status){
				 cout<<"Error closing file "<<HIR_filename<<". Status= "<<status<<endl;
				 return(status);
			 }
		 }
		if(galdef.verbose==-1003){// selectable debug
			strcpy(HIR_filename,configure.fits_directory);
			strcat(HIR_filename,"HIR_healpix.fits");// |b|>39.75  S.Digel
			SkymapToFits(galaxy.hpHIR, HIR_filename, "unit", "type");
		}
		if(galdef.verbose==-303)// selectable debug
        {
          cout<<"read_HIR: galaxy.HIR:"<<endl;
		      galaxy.hpHIR.print(cout);
        }
	 }else{
		 */

   int NAXIS,NAXIS1,NAXIS2,NAXIS3,NAXIS4;
   float CRVAL1,CRVAL2,CRVAL3;
   float CDELT1,CDELT2,CDELT3;
   char comment[100];

   Distribution HIR_input; // HI rings as read from FITS file
	 Skymap<double> nusehp;
   Distribution nuse;      // number of cells used in rebinned map

   //   cout<<"galaxy.n_long,n_lat  "<<galaxy.n_long<<" "<<galaxy.n_lat<<endl;

   strcpy(HIR_filename,configure.fits_directory);
   // strcat(HIR_filename,"COMPASS.HIR.M0010051");// 6 rings Ro=10 kpc
   if(galdef.HI_survey==8)                                                //AWS20050913
      strcat(HIR_filename,"COMPASS.HIR.M0010058");// 8 rings Digel Ro=8.5 kpc
   if(galdef.HI_survey==9)
//    strcat(HIR_filename,"rbands_hi2_32.fits");// 9 rings Digel Aug  2005  AWS20050826
//    strcat(HIR_filename,"rbands_hi5.fits"   );// 9 rings Digel Jan  2006  AWS20060112
      strcat(HIR_filename,"rbands_hi7.fits"   );// 9 rings Digel Jan  2006  AWS20060116

   cout<<" reading HIR from "<<HIR_filename<<endl;

   if( fits_open_file(&fptr,HIR_filename,READONLY,&status) ) cout<<"read HIR open status= "<<status<<endl;

   if( fits_read_key(fptr,TINT,"NAXIS" ,&NAXIS ,comment,&status) ) cout<<"0read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS1",&NAXIS1,comment,&status) ) cout<<"1read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS2",&NAXIS2,comment,&status) ) cout<<"2read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS3",&NAXIS3,comment,&status) ) cout<<"3read HIR status= "<<status<<endl;
 

   if( fits_read_key(fptr,TFLOAT,"CRVAL1",&CRVAL1,comment,&status) ) cout<<"5read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CRVAL2",&CRVAL2,comment,&status) ) cout<<"6read HIR status= "<<status<<endl;

   if( fits_read_key(fptr,TFLOAT,"CDELT1",&CDELT1,comment,&status) ) cout<<"8read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CDELT2",&CDELT2,comment,&status) ) cout<<"9read HIR status= "<<status<<endl;


   cout<<" NAXIS = "<<NAXIS <<endl;
   cout<<" NAXIS1,2,3   = "<<NAXIS1<<" "<<NAXIS2<<" "<<NAXIS3<<endl;
   cout<<" CRVAL1,2   = "<<CRVAL1<<" "<<CRVAL2<<endl;
   cout<<" CDELT1,2   = "<<CDELT1<<" "<<CDELT2<<endl;

   long nelements=NAXIS1*NAXIS2*NAXIS3, felement=1;
   float *HIR_in=new float[nelements];
   float nulval=1e-15;
   int anynul;

   if( fits_read_img(fptr,TFLOAT,felement,nelements,&nulval,HIR_in,&anynul,&status) )
      cout<<"#read HIR status= "<<status<<endl;

   // for(int i=0; i<nelements; i++) cout<<HIR_in[i]<<" ";

   cout<<"generating galaxy.HIR:"<<endl;


         int i_long,i_lat,i_Ring;
         int i_long_in,i_lat_in;
         int n_long_in=NAXIS1;
         int n_lat_in =NAXIS2;
         int n_Ring=NAXIS3;

         
      
				 //The gas rings store the values of the pixel boundaries while we want
				 //the average value for the pixel.  Therefore the -1 in latitude.
         HIR_input.init(n_long_in,n_lat_in-1,n_Ring,1);

         cout<<" HIR_input initialized"<<endl;



         for  ( i_Ring=0;  i_Ring<n_Ring; i_Ring++){
					 const int ind1 = i_Ring*n_lat_in*n_long_in;
	 
	  for ( i_lat_in =0;  i_lat_in <n_lat_in-1;  i_lat_in++ ){
			 const int i_lat_max = i_lat_in+1;
			const int ind2l = ind1 + i_lat_in*n_long_in;
			const int ind2u = ind1 + i_lat_max*n_long_in;
             
	   for( i_long_in=0;  i_long_in<n_long_in; i_long_in++){
			 const int i_long_max = i_long_in == (n_long_in -2) ? 0 : i_long_in +1;
            

                  HIR_input.d3[i_long_in][i_lat_in][i_Ring].s[0]= 
										(HIR_in[ind2l+i_long_in]+HIR_in[ind2u+i_long_in]+
										HIR_in[ind2l+i_long_max]+HIR_in[ind2u+i_long_max])/4.;

          
            }  //  i_long_in
           }   //  i_lat_in
	 }     //  i_Ring
   

         HIR_input*= 1.0e20;   // applies to HIR.10058, not HIR.10051   

	 if(galdef.verbose== -201) // selectable debug
      {
      
                HIR_input.print();
      }
  


	 if (galdef.skymap_format == 3){
	   galaxy.hpHIR.Resize(galdef.healpix_order, n_Ring, NEST);
		 galaxy.hpHIR.filllbCARarray(  HIR_input );
	 }else{

double l,b;

// create column density map corresponding to required gamma skymaps
// by rebinning input map
// cout<<"galay.HIR.init( "<<galaxy.n_long<<" "<<galaxy.n_lat<<" "<<n_Ring<<",1)"<<endl;

      galaxy.HIR.init(galaxy.n_long,galaxy.n_lat,n_Ring,1);

      nuse.init      (galaxy.n_long,galaxy.n_lat,n_Ring,1);


         for  (i_Ring=0;     i_Ring   <n_Ring; i_Ring++)   {
       	  for (i_lat_in =0;  i_lat_in <n_lat_in-1;  i_lat_in++ ){
           for(i_long_in=0;  i_long_in<n_long_in; i_long_in++){
            
						 //We have already moved from edge representation to average
						 //representation so we need to have the l and b values in the
						 //center of the pixel.
                  l=CRVAL1+CDELT1/2.+(i_long_in)*CDELT1;
                  b=CRVAL2+CDELT2/2.+(i_lat_in )*CDELT2;
									//Find the pixel that contains the l and b values of the
									//center of the pixel.  Due to the nature of int conversion
									//flooring the value, we use the lower left corner of the
									//pixel rather than the center.
                  i_long= (int) ((l-galaxy.long_min+galaxy.d_long/2.)/galaxy.d_long); // IMOS20010816
                  i_lat = (int) ((b-galaxy.lat_min +galaxy.d_lat/2. )/galaxy.d_lat);  // IMOS20010816

											if(i_long>=0 && i_long<galaxy.n_long&&i_lat >=0 && i_lat<galaxy.n_lat){

											galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]+= HIR_input.d3[i_long_in][i_lat_in][i_Ring].s[0];
											nuse.d3      [i_long][i_lat][i_Ring].s[0]+=1;
											}


		  //       cout<<"l b i_long_in i_lat_in i_long i_lat i_Ring "<<l<<" "<< b<<" "<< i_long_in<<" "<< i_lat_in<<" "<< i_long<<" "<< i_lat<<" "<< i_Ring<<endl;
                  

          
            }  //  i_long
           }   //  i_lat
	 }     //  i_Ring

}//HEALPIX OR NOT

	 

	/////////////////////////////////// |b|>39.75

  if(galdef.HI_survey==8)// only for 8 ring survey; 9 ring survey is all-sky  AWS20050914
  {

   strcpy(HIR_filename,configure.fits_directory);
   //   strcat(HIR_filename,"COMPASS.HIR.M0010049");// |b|>24.75, *0.5
   strcat(HIR_filename,"COMPASS.HIR.M0010055");// |b|>39.75  S.Digel
     
   cout<<" reading HIR from "<<HIR_filename<<endl;

   if( fits_open_file(&fptr,HIR_filename,READONLY,&status) ) cout<<"read HIR open status= "<<status<<endl;

   if( fits_read_key(fptr,TINT,"NAXIS" ,&NAXIS ,comment,&status) ) cout<<"0read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS1",&NAXIS1,comment,&status) ) cout<<"1read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS2",&NAXIS2,comment,&status) ) cout<<"2read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS3",&NAXIS3,comment,&status) ) cout<<"3read HIR status= "<<status<<endl;
 

   if( fits_read_key(fptr,TFLOAT,"CRVAL1",&CRVAL1,comment,&status) ) cout<<"5read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CRVAL2",&CRVAL2,comment,&status) ) cout<<"6read HIR status= "<<status<<endl;
;
   if( fits_read_key(fptr,TFLOAT,"CDELT1",&CDELT1,comment,&status) ) cout<<"8read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CDELT2",&CDELT2,comment,&status) ) cout<<"9read HIR status= "<<status<<endl;


   cout<<" NAXIS = "<<NAXIS <<endl;
   cout<<" NAXIS1,2,3   = "<<NAXIS1<<" "<<NAXIS2<<" "<<NAXIS3<<endl;
   cout<<" CRVAL1,2     = "<<CRVAL1<<" "<<CRVAL2<<endl;
   cout<<" CDELT1,2     = "<<CDELT1<<" "<<CDELT2<<endl;

        nelements=NAXIS1*NAXIS2*NAXIS3;
				delete[] HIR_in;          //Gulli20070810
          HIR_in=new float[nelements];
    

   if( fits_read_img(fptr,TFLOAT,felement,nelements,&nulval,HIR_in,&anynul,&status) )
      cout<<"#read HIR status= "<<status<<endl;

         n_long_in=NAXIS1;
         n_lat_in =NAXIS2;
      
				 HIR_input.delete_array();          //Gulli20070810
         HIR_input.init(n_long_in,n_lat_in-1,n_Ring,1);

         cout<<" HIR_input initialized for high latitudes"<<endl;



 
	 //         i_Ring=2            ;// local ring (0-4-8-10-12-15-30)//HIR-10051

	                       //                0   1   2   3V    4    5    6    7
          i_Ring=3            ;// local ring (1.5-3.5-5.5-7.5 - 9.5-11.5-13.5-15.5-50) HIR-10055
					 const int ind1 = i_Ring*n_lat_in*n_long_in;
	 
	  for ( i_lat_in =0;  i_lat_in <n_lat_in-1;  i_lat_in++ ){
			 const int i_lat_max = i_lat_in+1;
			const int ind2l = ind1 + i_lat_in*n_long_in;
			const int ind2u = ind1 + i_lat_max*n_long_in;
             
	   for( i_long_in=0;  i_long_in<n_long_in; i_long_in++){
			 const int i_long_max = i_long_in == (n_long_in -2) ? 0 : i_long_in +1;
            

                  HIR_input.d3[i_long_in][i_lat_in][i_Ring].s[0]= 
										(HIR_in[ind2l+i_long_in]+HIR_in[ind2u+i_long_in]+
										HIR_in[ind2l+i_long_max]+HIR_in[ind2u+i_long_max])/4.;


          
            }  //  i_long_in
           }   //  i_lat_in
 
          HIR_input*= 1.0e20;   // applies to HIR.10055, not HIR.10049   

	  if(galdef.verbose== -202) // selectable debug
      {
      
                HIR_input.print();
      }

	 if (galdef.skymap_format == 3){
		 Skymap<double> hpHIR_in(galaxy.hpHIR);
		 hpHIR_in.filllbCARarray(  HIR_input );
		 galaxy.hpHIR += hpHIR_in;
	 }else{
    

       	  for (i_lat_in =0;  i_lat_in <n_lat_in-1;  i_lat_in++ ){
           for(i_long_in=0;  i_long_in<n_long_in; i_long_in++){
            
						 //We have already moved from edge representation to average
						 //representation so we need to have the l and b values in the
						 //center of the pixel.
                  double l=CRVAL1+CDELT1/2.+(i_long_in)*CDELT1;
                  double b=CRVAL2+CDELT2/2.+(i_lat_in )*CDELT2;
                  if(b<-39.75 || b>+39.75){
									//Find the pixel that contains the l and b values of the
									//center of the pixel.  Due to the nature of int conversion
									//flooring the value, we use the lower left corner of the
									//pixel rather than the center.
                  i_long= (int) ((l-galaxy.long_min+galaxy.d_long/2.)/galaxy.d_long); // IMOS20010816
                  i_lat = (int) ((b-galaxy.lat_min +galaxy.d_lat/2. )/galaxy.d_lat);  // IMOS20010816



                  if(i_long>=0 && i_long<galaxy.n_long&&i_lat >=0 && i_lat<galaxy.n_lat){

	                                                            
                  galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]+= HIR_input.d3[i_long_in][i_lat_in][i_Ring].s[0]    ;
                  nuse.d3      [i_long][i_lat][i_Ring].s[0]+=1;
									}
									}

		  //	         cout<<"l b i_long_in i_lat_in i_long i_lat i_Ring "<<l<<" "<< b<<" "<< i_long_in<<" "<< i_lat_in<<" "<< i_long<<" "<< i_lat<<" "<< i_Ring<<endl;
                  

          
            }  //  i_long
           }   //  i_lat
	 } //HEALPIX OR NOT
       

  }// if (HI_survey==8)

	//////////////////////////////////////////////////

	 

	 if (galdef.skymap_format == 3){
		if(galdef.verbose==-1003){// selectable debug
			strcpy(HIR_filename,configure.fits_directory);
			strcat(HIR_filename,"HIR_healpix.fits");// |b|>39.75  S.Digel
			SkymapToFits(galaxy.hpHIR, HIR_filename, "unit", "type");
		}
		if(galdef.verbose==-303)// selectable debug
        {
          cout<<"read_HIR: galaxy.HIR:"<<endl;
		      galaxy.hpHIR.print(cout);
        }
	 }else{
	 // normalize by number of cells used in rebinning
        for  (i_Ring=0;  i_Ring<       n_Ring; i_Ring++){
	 for (i_lat =0;  i_lat <galaxy.n_lat;  i_lat++ ){
          for(i_long=0;  i_long<galaxy.n_long; i_long++){

            if(nuse.d3      [i_long][i_lat][i_Ring].s[0]>0)
               galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]/=nuse.d3      [i_long][i_lat][i_Ring].s[0];
          
            }  //  i_long
           }   //  i_lat
	}//i_Ring
       nuse.delete_array();
       if(galdef.verbose== -203)// selectable debug
      {
        cout<<"read_HIR: galaxy.HIR:"<<endl;
                galaxy.HIR.print();
      }
	 } //NOT HEALPIX


       delete[] HIR_in;
       HIR_input.delete_array();  
//}


	//////////////////////////////////////////////////

   cout<<" <<<< read_HIR    "<<endl;
   return status;
}
