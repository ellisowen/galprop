using namespace std;
#include<iostream>


      void inter(double* XX,double* FF,int NN, double X0, double &F1)
{
	/*
      SUBROUTINE DINTER(XX,FF,NN,X0,F1)

C********************************************************************
C* INTERPOLATION BY CUBIC POLINOM, REAL*8 version                   *
C* XX,FF  ARE THE ARGUMENT AND FUNCTION ARRAYS (j=1..NN);           *
C* X0  IS THE ARGUMENT AT WHICH THE FUNCTION IS TO BE EVALUATED;    *
C* F1  ON EXIT, IS THE EVALUATED VALUE OF THE FUNCTION.             *
C********************************************************************
      IMPLICIT real*8  (A-H,O-Z)
      IMPLICIT integer (I-N)



      DIMENSION XX(NN),FF(NN)
keep the fortran indexing in this version for use with bremsstrahlung.cc

	*/

  int debug=0;
  if(debug==1)cout<<"DINTER"<<endl;

      int IN=0;
      int N1=NN-1;

      int M;
      for(    M=1;M<=N1;M++)if(X0>= XX[M] &&  X0 <= XX[M+1]) break; 
           


      IN=M;
      if(IN ==  1 || X0 <  XX[1] ) IN=2;
      if(IN == N1 || X0<   XX[NN]) IN=NN-2;

       double X1=XX[IN-1];
       double X2=XX[IN];
       double X3=XX[IN+1];
       double X4=XX[IN+2];
       double Y1=FF[IN-1];
       double Y2=FF[IN];
       double Y3=FF[IN+1];
       double Y4=FF[IN+2];

       double X01=X0-X1  ;          
       double X02=X0-X2  ;    
       double X03=X0-X3  ;    
       double X04=X0-X4  ;    
       double X12=X1-X2 ;     
       double X13=X1-X3   ;   
       double X14=X1-X4   ;   
       double X23=X2-X3   ;   
       double X24=X2-X4   ;   
       double X34=X3-X4;

      X1= X02*X03*X04/X12/X13/X14;
      X2=-X01*X03*X04/X12/X23/X24;
      X3= X01*X02*X04/X13/X23/X34;
      X4=-X01*X02*X03/X14/X24/X34;

      F1= Y1*X1+Y2*X2+Y3*X3+Y4*X4;

      return    ;
      
     };
