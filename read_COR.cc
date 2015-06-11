
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * read_COR.cc *                                galprop package * 5/22/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"
#include"fitsio.h"
#include"Skymap.h"
#include <vector>
#include <algorithm>
#include <string.h>
#include"SkymapFitsio.h"
#include<cmath>

int Galprop::read_COR()
{
   cout<<" >>>> read_COR    "<<endl;
   int status=0;
   fitsfile *fptr;
   char COR_filename[200];

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
	   galaxy.hpCOR.Resize(galdef.healpix_order, n_Ring, NEST);
		 //Number of pixels in read in map that is equivalent to one pixel in the
		 //used map
		 int nside = 1<<order;
		 long npix = 12*nside*nside;
		 const int pixelRatio = npix/galaxy.hpCOR.Npix();
		 const int invPixelRatio = galaxy.hpCOR.Npix()/npix;
		 //The fits directory
		 char COR_file[200];
		 for (int i_Ring=0; i_Ring < n_Ring; ++i_Ring){
			 sprintf(COR_file, "rbands_co4_%d_512.fits", i_Ring);
			 strcpy(COR_filename,configure.fits_directory);
			 strcat(COR_filename,COR_file);
			 if( fits_open_file(&fptr,COR_filename,READONLY,&status) ) {
				 cout<<"read COR open status= "<<status<<endl;
				 cout<<"COR filename= "<<COR_filename<<endl;
				 return(status);
			 }
			 //Move to the correct HDU number 2
			 if( fits_movabs_hdu(fptr,2,&hdutype,&status) ) {
				 cout<<"move to COR binary table status= "<<status<<endl;
				 return status;
			 }
			 //How we read the map depends on which map has the finer grid.
			 //pixelRatio will be zero if hpCOR has finer grid
			 if (pixelRatio > 1){
				 double *tmparray = new double[pixelRatio];
				 //Loop through the hpCOR map and add it one pixel at a time
				 for( int p=0; p<galaxy.hpCOR.Npix(); ++p ){
					 long frow = p*pixelRatio/nrow + 1;
					 long felement = p*pixelRatio%nrow + 1;
					 fits_read_col(fptr,TDOUBLE,1,frow,felement,pixelRatio,&nulval,tmparray,&anynul,&status);
					 if (status){
						 cout<<"read COR read table failed, status= "<<status<<endl;
						 cout<<"p= "<<p<<", pixelRatio= "<<pixelRatio<<", felement= "<<felement<<endl;
						 return(status);
					 }
					 //Calculate the average and insert into skymap
					 double sum = 0;
					 for( int i=0; i<pixelRatio; ++i){
						 sum += tmparray[i];
					 }
					 galaxy.hpCOR[p][i_Ring] = sum/double(pixelRatio);
				 }
				 delete[] tmparray;
			 }else{
				 double *tmparray = new double[npix];
				 //Read the whole skymap
				 fits_read_col(fptr,TDOUBLE,1,1L,1L,npix,&nulval,tmparray,&anynul,&status);
				 if (status){
					 cout<<"read COR read table failed, status= "<<status<<endl;
					 return(status);
				 }
				 //Loop through the array and insert into 
				 for(int i=0; i<npix; ++i){
					 for(int p=i*invPixelRatio; p<(i+1)*invPixelRatio; ++p){
						 galaxy.hpCOR[p][i_Ring] = tmparray[i];
					 }
				 }
				 delete[] tmparray;
			 }
			 fits_close_file(fptr, &status);
			 if(status){
				 cout<<"Error closing file "<<COR_filename<<". Status= "<<status<<endl;
				 return(status);
			 }
		 }
		if(galdef.verbose==-1003){// selectable debug
			strcpy(COR_filename,configure.fits_directory);
			strcat(COR_filename,"COR_healpix.fits");// |b|>39.75  S.Digel
			SkymapToFits(galaxy.hpCOR, COR_filename, "unit", "type");
		}
		if(galdef.verbose==-303)// selectable debug
        {
          cout<<"read_COR: galaxy.COR:"<<endl;
		      galaxy.hpCOR.print(cout);
        }
	 }else{
		 */

   int NAXIS,NAXIS1,NAXIS2,NAXIS3,NAXIS4;
   float CRVAL1,CRVAL2,CRVAL3;
   float CDELT1,CDELT2,CDELT3;
   char comment[100];

   Distribution COR_input; // CO rings as read from FITS file
   Distribution nuse;      // number of cells used in rebinned map


   strcpy(COR_filename,configure.fits_directory);
   //   strcat(COR_filename,"COMPASS.COR.M0010453");// 6 rings Ro=10 kpc
   if(galdef.CO_survey==8)                     //AWS20050913
   strcat(COR_filename,"COMPASS.COR.M0010465");// 8 rings Ro=8.5 kpc Digel
   if(galdef.CO_survey==9)                     //AWS20050913
// strcat(COR_filename,"COMPASS.COR.M0010465_9rings");// 8 rings Ro=8.5 kpc Digel, modified to 9 ring format
// strcat(COR_filename,"rbands_co2.fits"            );// 9 rings Ro=8.5 kpc Digel 13 Jan 2006  //AWS20060113
// strcat(COR_filename,"rbands_co3.fits"            );// 9 rings Ro=8.5 kpc Digel 16 Jan 2006  //AWS20060116
   strcat(COR_filename,"rbands_co4.fits"            );// 9 rings Ro=8.5 kpc Digel 18 Jan 2006  //AWS20060118

   cout<<" reading COR from "<<COR_filename<<endl;

   if( fits_open_file(&fptr,COR_filename,READONLY,&status) ) cout<<"read COR open status= "<<status<<endl;

   if( fits_read_key(fptr,TINT,"NAXIS" ,&NAXIS ,comment,&status) ) cout<<"0read COR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS1",&NAXIS1,comment,&status) ) cout<<"1read COR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS2",&NAXIS2,comment,&status) ) cout<<"2read COR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS3",&NAXIS3,comment,&status) ) cout<<"3read COR status= "<<status<<endl;


   if( fits_read_key(fptr,TFLOAT,"CRVAL1",&CRVAL1,comment,&status) ) cout<<"5read COR status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CRVAL2",&CRVAL2,comment,&status) ) cout<<"6read COR status= "<<status<<endl;
 
   if( fits_read_key(fptr,TFLOAT,"CDELT1",&CDELT1,comment,&status) ) cout<<"8read COR status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CDELT2",&CDELT2,comment,&status) ) cout<<"9read COR status= "<<status<<endl;
 

   cout<<" NAXIS = "<<NAXIS <<endl;
   cout<<" NAXIS1,2,3 = "<<NAXIS1<<" "<<NAXIS2<<" "<<NAXIS3<<endl;
   cout<<" CRVAL1,2   = "<<CRVAL1<<" "<<CRVAL2<<endl;
   cout<<" CDELT1,2   = "<<CDELT1<<" "<<CDELT2<<endl;

   long nelements=NAXIS1*NAXIS2*NAXIS3, felement=1;
   float *COR_in=new float[nelements];
   float nulval=1e-15;
   int anynul;

   if( fits_read_img(fptr,TFLOAT,felement,nelements,&nulval,COR_in,&anynul,&status) )
      cout<<"#read COR status= "<<status<<endl;

   // for(int i=0; i<nelements; i++) cout<<COR_in[i]<<" ";

   cout<<"generating galaxy.COR:"<<endl;


         int i_long,i_lat,i_Ring;
         int i_long_in,i_lat_in;
         int n_long_in=NAXIS1;
         int n_lat_in =NAXIS2;
         int n_Ring=NAXIS3;

         
      
				 //The gas rings store the values of the pixel boundaries while we want
				 //the average value for the pixel.  Therefore the -1 in latitude.
         COR_input.init(n_long_in,n_lat_in-1,n_Ring,1);

         cout<<" COR_input initialized"<<endl;


         for  ( i_Ring=0;  i_Ring<n_Ring; i_Ring++){
					 const int ind1 = i_Ring*n_lat_in*n_long_in;
	 
	  for ( i_lat_in =0;  i_lat_in <n_lat_in-1;  i_lat_in++ ){
			 const int i_lat_max = i_lat_in+1;
			const int ind2l = ind1 + i_lat_in*n_long_in;
			const int ind2u = ind1 + i_lat_max*n_long_in;
             
	   for( i_long_in=0;  i_long_in<n_long_in; i_long_in++){
			 const int i_long_max = i_long_in == (n_long_in -2) ? 0 : i_long_in +1;
            

                  COR_input.d3[i_long_in][i_lat_in][i_Ring].s[0]= 
										(COR_in[ind2l+i_long_in]+COR_in[ind2u+i_long_in]+
										COR_in[ind2l+i_long_max]+COR_in[ind2u+i_long_max])/4.;
	 
            }  //  i_long_in
           }   //  i_lat_in
	 }     //  i_Ring
   



	 if(galdef.verbose== -301 )//selectable debug
      {
      
                COR_input.print();
      }
  

	 if (galdef.skymap_format == 3){
	   galaxy.hpCOR.Resize(galdef.healpix_order, n_Ring, NEST);
		 galaxy.hpCOR.filllbCARarray(  COR_input );
		if(galdef.verbose==-1003){// selectable debug
			strcpy(COR_filename,configure.fits_directory);
			strcat(COR_filename,"COR_healpix.fits");
			SkymapToFits(galaxy.hpCOR, COR_filename, "unit", "type");
		}
		if(galdef.verbose==-303)// selectable debug
        {
          cout<<"read_COR: galaxy.COR:"<<endl;
		      galaxy.hpCOR.print(cout);
        }
	 }else{

double l,b;

// create column density map corresponding to required gamma skymaps
// by rebinning input map
// An extra ring is added for compatibility with HI rings (6 instead of the input 5) [COR-100453 only]
// and this CO ring is left at zero.

      galaxy.COR.init(galaxy.n_long,galaxy.n_lat,n_Ring,1);

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

                  galaxy.COR.d3[i_long][i_lat][i_Ring].s[0]+= COR_input.d3[i_long_in][i_lat_in][i_Ring].s[0];
                  nuse.d3      [i_long][i_lat][i_Ring].s[0]+=1;


		  //       cout<<"l b i_long_in i_lat_in i_long i_lat i_Ring "<<l<<" "<< b<<" "<< i_long_in<<" "<< i_lat_in<<" "<< i_long<<" "<< i_lat<<" "<< i_Ring<<endl;
		  }
                  

          
            }  //  i_long
           }   //  i_lat
	 }     //  i_Ring

	 // normalize by number of cells used in rebinning
        for  (i_Ring=0;  i_Ring<       n_Ring; i_Ring++){
	 for (i_lat =0;  i_lat <galaxy.n_lat;  i_lat++ ){
          for(i_long=0;  i_long<galaxy.n_long; i_long++){

            if(nuse.d3      [i_long][i_lat][i_Ring].s[0]>0)
               galaxy.COR.d3[i_long][i_lat][i_Ring].s[0]/=nuse.d3      [i_long][i_lat][i_Ring].s[0];
          
            }  //  i_long
           }   //  i_lat
	 }     //  i_Ring
       nuse.delete_array();
	if(galdef.verbose==-303)// selectable debug
      {
        cout<<"read_COR: galaxy.COR:"<<endl;
                galaxy.COR.print();
      }
   cout<<" <<<< read_COR    "<<endl;
	 }



       delete[] COR_in;
       COR_input.delete_array();  



//}

   return status;
}

