//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * GasCubeData.h *                                   galprop package * 08/16/2001
// * Author Eric Carlson 10/26/2014
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
#ifndef _GasCubeData_h_
#define _GasCubeData_h_

using namespace std;

struct GasCubeData{
    float *Data;
    double NAXIS1, NAXIS2, NAXIS3;
    double CRVAL1, CRVAL2, CRVAL3;
    double CDELT1, CDELT2, CDELT3;
    double total; 
};

#endif