/**************************************************************************


     DRIVE DAS SUBROTINAS

     double viscosity (double t, double ro)

 *************************************************************************/


#include <mex.h>

#include <math.h>
#include <stdio.h>

int  lll[4], alfaz[4];
unsigned betaz[4];
double  tttc, pppc, vvvc, rrroc, mmmolwt, ggg[4], bbb[10][7], aaa[5][2];


#define ro_inf_init 0.001
#define ro_sup_init 1.50
#define perc_er 1.0e-6
#define it_max 100



double viscosity (double t, double ro)
{
   double w, tau, visc;
   int i, j;
   static int readparameters=0;
   void parameters(void);

   if (readparameters==0)
   {
      readparameters=1;
      parameters();
   }

   visc=0.;
   w=ro/0.468;
   tau=t/304.2;

   for (i=1 ; i<5 ; i++)
   {
     for (j=0 ; j<2 ; j++)
	 visc+=aaa[i][j]*pow(w,i)/pow(tau,j);
   }

   visc=exp(visc);
   visc*=sqrt(tau)*(2722.46461-1663.46068/tau+466.920556/tau/tau)*1e-7;

   return visc;
}




/* *************************************************************************
   ************************************************************************ */


void parameters()
{
  // CARBON DIOXIDE CRITICAL CONSTANTES
     tttc=304.1;       // K
     pppc=72.8349174;  // atm
     vvvc=93.9;        // cm3/mol
     rrroc=0.4686901;  // g/cm3
     mmmolwt=44.010;   // g/mol

  // PARAMETERS OF THE ALTUNIN-SAKHABETINOV
     aaa[1][0]= 0.248566120;
     aaa[1][1]= 0.004894942;
     aaa[2][0]=-0.373300660;
     aaa[2][1]= 1.227534880;
     aaa[3][0]= 0.363854523;
     aaa[3][1]=-0.774229021;
     aaa[4][0]=-0.0639070755;
     aaa[4][1]=0.1425070490;

  // PARAMETERS OF THE AG-HGK EQUATION
     lll[1]=0;
     lll[2]=2;
     lll[3]=0;
     alfaz[1]=34;
     alfaz[2]=40;
     alfaz[3]=30;
     betaz[1]=20000;
     betaz[2]=20000;
     betaz[3]=40000;
     ggg[1]=-7.53e-4;
     ggg[2]=-5.73e-3;
     ggg[3]= 1.84e-4;
     bbb[0][0]=-0.7255896770;
     bbb[0][1]=-1.6698566330;
     bbb[0][2]= 0.4191613578;
     bbb[0][3]= 1.1540585470;
     bbb[0][4]= 1.1450275820;
     bbb[0][5]= 1.1488455130;
     bbb[0][6]= 0.7069388840;
     bbb[1][0]= 0.4481451002;
     bbb[1][1]= 1.2690839330;
     bbb[1][2]= 6.0578119110;
     bbb[1][3]= 15.859789780;
     bbb[1][4]= 20.218370270;
     bbb[1][5]= 9.1900771440;
     bbb[2][0]=-0.1743673384;
     bbb[2][1]=-1.9544044470;
     bbb[2][2]=-5.6151979650;
     bbb[2][3]=-6.9768169150;
     bbb[2][4]=-0.5761694929;
     bbb[2][5]= 3.0072849370;
     bbb[3][0]=-4.2438160930e-4;
     bbb[3][1]=-1.7884558440;
     bbb[3][2]=-11.346293670;
     bbb[3][3]=-29.104035620;
     bbb[3][4]=-30.026639370;
     bbb[3][5]=-8.3612823860;
     bbb[4][0]= 0.2668130548;
     bbb[4][1]= 2.7185742230;
     bbb[4][2]= 9.4622888160;
     bbb[4][3]= 10.603173790;
     bbb[4][4]= 0.1567993789;
     bbb[4][5]=-2.7232168500;
     bbb[5][0]= 0.07340283381;
     bbb[5][1]= 1.1547892190;
     bbb[5][2]= 7.4509888050;
     bbb[5][3]= 16.001430470;
     bbb[5][4]= 10.971048690;
     bbb[6][0]=-0.1756082074;
     bbb[6][1]=-2.1141845860;
     bbb[6][2]=-6.1447687020;
     bbb[6][3]=-4.6675661180;
     bbb[7][0]= 8.8442710160e-3;
     bbb[7][1]= 0.0148894556;
     bbb[7][2]=-1.4450102070;
     bbb[7][3]=-1.9979431860;
     bbb[8][0]= 0.06107749242;
     bbb[8][1]= 0.6239980516;
     bbb[8][2]= 1.1940662950;
     bbb[9][0]=-0.01994277669;
     bbb[9][1]=-0.1666138543;
     bbb[9][2]= 5.9238882890e-3;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{           
  double viscosity_calculada;
  double *auxPtr;
  double T;
  double rho;
  
  T = mxGetScalar(prhs[0]); // primeira possi��o
  rho = mxGetScalar(prhs[1]); // segunda possi��o

  viscosity_calculada = (double)viscosity(T,rho); 
  
  /* Create a 1-by-1 matrix for the return argument. */
   plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
   auxPtr = (double*)mxGetPr(plhs[0]);
   *auxPtr = viscosity_calculada;
   
  return;
                  
}


