
#include "../paul.h"

void prim2cons( double * , double * , double , double , double * );
void cons2prim( double * , double * , double , double , double * );

double get_pre( double * , double * );
double get_temp( double * , double * );

static int NL = 0;
static double * rr  = NULL;
static double * rho = NULL;
static double * Pp  = NULL;
static double * vr  = NULL;
static double * Om  = NULL;
//static double * Ss  = NULL;

int countlines(char * filename){
   FILE *pFile = fopen(filename, "r");
   int lines=0;
   char c;
   while ((c = fgetc(pFile)) != EOF){
      if (c == '\n') ++lines;
   }
   fclose(pFile);
   return(lines);
}

int getTable( void ){
   int nL = countlines("Initial/initial.dat");
   int error;
   rr  = (double *) malloc( nL*sizeof(double) );
   rho = (double *) malloc( nL*sizeof(double) );
   Pp  = (double *) malloc( nL*sizeof(double) );
   vr  = (double *) malloc( nL*sizeof(double) );
   Om  = (double *) malloc( nL*sizeof(double) );
   //Ss  = (double *) malloc( nL*sizeof(double) );
   FILE * pFile = fopen("Initial/initial.dat","r");
   int l;
   for( l=0 ; l<nL ; ++l ){
      //error = fscanf(pFile,"%le %le %le %le %le %le %le\n",&(rr[l]),&(rho[l]),&(Pp[l]),&(vr[l]),&(Om[l]),&(Ss[l]));
      error = fscanf(pFile,"%le %le %le %le %le \n",&(rr[l]),&(rho[l]),&(Pp[l]),&(vr[l]),&(Om[l]));
   }
   fclose(pFile);
   return(nL);
}

void setICparams( struct domain * theDomain ){
   NL = getTable();
}

void initial( double * prim , double r , double * T ){

   int l=0;
   while( rr[l] < r && l < NL-2 ) ++l;
   if( l==0 ) ++l;

   double rp = rr[l];
   double rm = rr[l-1];
   double drm = fabs(r-rm);
   double drp = fabs(rp-r);

   double rh    = (rho[l-1]*drp + rho[l]*drm)/(drp+drm);
   double P     = (Pp[l-1]*drp  + Pp[l]*drm )/(drp+drm);
   double V     = (vr[l-1]*drp  + vr[l]*drm )/(drp+drm);
   double X     = (Om[l-1]*drp  + Om[l]*drm )/(drp+drm);
   //double S     = (Ss[l-1]*drp  + Ss[l]*drm )/(drp+drm);

   /*if( l==NL-2 ){
      rh = rho[l]; 
      P  = Pp[l];
      V  = vr[l];
      X  = Om[l];
   }*/

   double x[NUM_I];
   x[0] = 0.0;          //p
   x[1] = X;            //He4
   x[2] = 0.5*(1.0-X);  //C12
   x[3] = 0.0;          //N13
   x[4] = 0.0;          //N14
   x[5] = 0.5*(1.0-X);  //O16
   x[6] = 0.0;          //F18
   x[7] = 0.0;          //Ne20
   x[8] = 0.0;          //Ne21
   x[9] = 0.0;          //Na22
   x[10] = 0.0;         //Na23
   x[11] = 0.0;         //Mg24
   x[12] = 0.0;         //Al27
   x[13] = 0.0;         //Si28
   x[14] = 0.0;         //P31
   x[15] = 0.0;         //S32
   x[16] = 0.0;         //Ar36
   x[17] = 0.0;         //Ca40
   x[18] = 0.0;         //Ti44
   x[19] = 0.0;         //Cr48
   x[20] = 0.0;         //Fe52
   x[21] = 0.0;         //Ni56
   
   prim[RHO] = rh;
   prim[PPP] = P;
   prim[VRR] = V;
   prim[AAA] = 0.0;
   for( int i=0 ; i<NUM_I ; i++ ){ prim[XXX+i] = x[i]; }

   //////////////////////// GET temperature ///////////////////////////
   double temp = 1e7; //a guess
   temp = get_temp( prim , &temp );
   *T = temp;

   //IGNITE CELLS//  
   double r_ignition = 4.5e8;//0.0;
   double width = 5e6;//2e6;
   double temp_peak = 2e9;//2e9;
   if( fabs(rr[l]-r_ignition)<2*width ){
      temp = temp+(temp_peak-temp)*pow(M_PI,-(rr[l]-r_ignition)*(rr[l]-r_ignition)/(.25*width*width));
      prim[PPP] = get_pre( prim , &temp );
      *T = temp;
   }
}

void initial_grid( int ii , double * rr_i , double * dr_i ){

}

void freeTable( void ){
   free(rr);
   free(rho);
   free(Pp);
   free(vr);
   free(Om);
   //free(Ss);
}
