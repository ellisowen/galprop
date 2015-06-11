
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * read_gas_cube.cc *                             galprop package * 10/24/2014 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
//written by Eric C. Carlson (ECC20141024)

using namespace std;
#include"galprop_classes.h"
#include"galprop_internal.h"
#include"fitsio.h"
#include"Skymap.h"
//#include <vector>
//#include <algorithm>
//#include"SkymapFitsio.h"
#include <cstring>
#include <string>
#include "Galdef.h"

//class GasCube

/* Read in the Specified data cube (gas density as a function of position)*/
/* Read in the Specified data cube (gas density as a function of position)
   Authored by Eric C. Carlson (ECC20141024) */
int GasCube::read_gas_cubes(char* type)
 {
    cout<<" >>>> read_gas_cubes "<<type<<endl;
    int status=0;
    fitsfile *fptr;
    int NAXIS,NAXIS1,NAXIS2,NAXIS3;
    float CRVAL1,CRVAL2,CRVAL3;
    float CDELT1,CDELT2,CDELT3;
    char comment[100];//, filename[200];

    const std::string fitsDirectory = configure.fFITSDataDirectory;
    std::string filename = fitsDirectory;

    Distribution HICube_input;
    Distribution COCube_input;

    if( !strcmp ("HICube",type) ) filename += galdef.HICube_filename;//strcat(filename,galdef.HIR_filename);
    if( !strcmp ("COCube",type) ) filename += galdef.COCube_filename;//strcat(filename,galdef.COR_filename);

    cout<<" reading "<<type<<" from "<<filename<<endl;

    if( fits_open_file(&fptr,filename.c_str(),READONLY,&status) ) {
       cout<<"read "<<type<<" open status= "<<status<<endl;

       if (strcmp ("HICube",type)) {
        throw std::runtime_error("Failed to open HI data cube \"" + filename + "\". Check HICube_filename in galdef.");
       }
       if (strcmp ("COCube",type)) {
        throw std::runtime_error("Failed to open CO data cube \"" + filename + "\".  Check COCube_filename in galdef.");
       }
       //if (strcmp(" ",galdef.COCube_filename) && strcmp (COCube",type))  cout << "Error: "<< type << " data cube not specified in galdef"<<endl;
       
    }

    if( fits_read_key(fptr,TINT,"NAXIS" ,&NAXIS ,comment,&status) ) cout<<"0read "<<type<<" status= "<<status<<endl;
    if( fits_read_key(fptr,TINT,"NAXIS1",&NAXIS1,comment,&status) ) cout<<"1read "<<type<<" status= "<<status<<endl;
    if( fits_read_key(fptr,TINT,"NAXIS2",&NAXIS2,comment,&status) ) cout<<"2read "<<type<<" status= "<<status<<endl;
    if( fits_read_key(fptr,TINT,"NAXIS3",&NAXIS3,comment,&status) ) cout<<"3read "<<type<<" status= "<<status<<endl;
    
    if( fits_read_key(fptr,TFLOAT,"CRVAL1",&CRVAL1,comment,&status) ) cout<<"4read "<<type<<" status= "<<status<<endl;
    if( fits_read_key(fptr,TFLOAT,"CRVAL2",&CRVAL2,comment,&status) ) cout<<"5read "<<type<<" status= "<<status<<endl;
    if( fits_read_key(fptr,TFLOAT,"CRVAL3",&CRVAL3,comment,&status) ) cout<<"6read "<<type<<" status= "<<status<<endl;
    
    if( fits_read_key(fptr,TFLOAT,"CDELT1",&CDELT1,comment,&status) ) cout<<"7read "<<type<<" status= "<<status<<endl;
    if( fits_read_key(fptr,TFLOAT,"CDELT2",&CDELT2,comment,&status) ) cout<<"8read "<<type<<" status= "<<status<<endl;
    if( fits_read_key(fptr,TFLOAT,"CDELT3",&CDELT3,comment,&status) ) cout<<"9read "<<type<<" status= "<<status<<endl;
    
    cout<<" NAXIS      = "<<NAXIS <<endl;
    cout<<" NAXIS1,2,3 = "<<NAXIS1<<" "<<NAXIS2<<" "<<NAXIS3<<endl;
    cout<<" CRVAL1,2,3   = "<<CRVAL1<<" "<<CRVAL2<<" "<<CRVAL3<<endl;
    cout<<" CDELT1,2,3   = "<<CDELT1<<" "<<CDELT2<<" "<<CDELT3<<endl;

    // X_CO = 2.3e20 for CO cube
    // Data scale for HI is 1e20 

    
    long nelements=NAXIS1*NAXIS2*NAXIS3, felement=1;
    float *image=new float[nelements];
    float nulval=1e-15;
    int anynul;

    // Read the image data 
    int HDU_TYPE;

    if( fits_movabs_hdu(fptr,1,&HDU_TYPE,&status) ) cout<<"10read "<<type<<" status= "<<status<<endl;
    if( fits_read_img(fptr,TFLOAT,felement,nelements,&nulval,image,&anynul,&status) )
    cout<<"#read "<<type<<" status= "<<status<<endl;
    
    return status ;
  }