
#include "../paul.h"

void prim2cons( double * , double * , double , double , double * );
void cons2prim( double * , double * , double , double , double * );

double get_pre( double * , double * );
double get_temp( double * , double * );
double get_eint( double * , double * );
double get_pre_from_etot( double * , double * );

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
   x[0] = X; //0.9*X;
   x[1] = 0.5*(1.0-X); //0.5*(1.0-0.9*X);
   x[2] = 0.5*(1.0-X); //0.5*(1.0-0.9*X);
   x[3] = 0.0; 
   x[4] = 0.0;
   x[5] = 0.0;
   x[6] = 0.0; 
   x[7] = 0.0;
   x[8] = 0.0;
   x[9] = 0.0;
   x[10] = 0.0;
   x[11] = 0.0;
   x[12] = 0.0;
   
   prim[RHO] = rh;
   prim[PPP] = P;
   prim[VRR] = V;
   prim[AAA] = 0.0;
   prim[XXX]      = x[0];
   prim[XXX+1]    = x[1];
   prim[XXX+2]    = x[2];
   prim[XXX+3]    = x[3];
   prim[XXX+4]    = x[4];
   prim[XXX+5]    = x[5];
   prim[XXX+6]    = x[6];
   prim[XXX+7]    = x[7];
   prim[XXX+8]    = x[8];
   prim[XXX+9]    = x[9];
   prim[XXX+10]   = x[10];
   prim[XXX+11]   = x[11];
   prim[XXX+12]   = x[12];

   //////////////////////// GET temperature ///////////////////////////
   double temp = 1e7; //a guess
   temp = get_temp( prim , &temp );
   *T = temp;

   //IGNITE CELLS//  
   /*
   double r_ignition = 4.5e8;
   double width = 2.5e6;     //one gaussian sigma
   double temp_peak = 2e9;
   if( fabs(rr[l]-r_ignition)<4*width ){
      temp = temp+(temp_peak-temp)*pow(M_PI,-(rr[l]-r_ignition)*(rr[l]-r_ignition)/(2*width*width));
      prim[PPP] = get_pre( prim , &temp );
      *T = temp;
   }
   */
   /////////////////////////////////////////////////
   /*
   double r_ignition = 5.1e8;
   double pert_temp_factor = 30.0;
   double pert_rad_factor = 0.5 * 2.5e6;
   double r1 = fabs(rr[l]-r_ignition) / pert_rad_factor;
   double X_he = prim[XXX];
   temp = temp * (1.0 + X_he * pert_temp_factor * 0.150 * (1.0 + tanh(2.0 - r1)));
   prim[PPP] = get_pre( prim , &temp );
   *T = temp;
   */
   /////////////////////////////////////////////////
   
   int q;
   double etot = 0.0, ebind1 = 0.0, ebind2 = 0.0;
   etot = get_eint( prim , &temp );
   for( q=0 ; q<NUM_I ; ++q ){ ebind1 += EBIND[q]*prim[q+XXX]; }

   double r_ignition = 5e8;
   double width = 1.25e6;//2e6;
   double flip = 0.2;//20% He4 to C12
   double frac_flip = flip*pow(M_PI,-(rr[l]-r_ignition)*(rr[l]-r_ignition)/(2*width*width));
   if( fabs(rr[l]-r_ignition)<4*width ){
      //Flip to nickel
      //prim[XXX+12] += prim[XXX]*frac_flip;
      //prim[XXX+12] += prim[XXX+1]*frac_flip;
      //prim[XXX+12] += prim[XXX+2]*frac_flip;

      //prim[XXX] -= prim[XXX]*frac_flip;
      //prim[XXX+1] -= prim[XXX+1]*frac_flip;
      //prim[XXX+2] -= prim[XXX+2]*frac_flip;
      
      //Flip to silicon
      //prim[XXX+5] += prim[XXX]*frac_flip;
      //prim[XXX+5] += prim[XXX+1]*frac_flip;
      //prim[XXX+5] += prim[XXX+2]*frac_flip;

      //prim[XXX] -= prim[XXX]*frac_flip;
      //prim[XXX+1] -= prim[XXX+1]*frac_flip;
      //prim[XXX+2] -= prim[XXX+2]*frac_flip;

      //Flip to sulfur 
      //prim[XXX+6] += prim[XXX]*frac_flip;
      //prim[XXX+6] += prim[XXX+1]*frac_flip;
      //prim[XXX+6] += prim[XXX+2]*frac_flip;

      //prim[XXX] -= prim[XXX]*frac_flip;
      //prim[XXX+1] -= prim[XXX+1]*frac_flip;
      //prim[XXX+2] -= prim[XXX+2]*frac_flip;

      //Flip to neon 
      prim[XXX+3] += prim[XXX]*frac_flip;
      prim[XXX+3] += prim[XXX+1]*frac_flip;
      prim[XXX+3] += prim[XXX+2]*frac_flip;

      prim[XXX] -= prim[XXX]*frac_flip;
      prim[XXX+1] -= prim[XXX+1]*frac_flip;
      prim[XXX+2] -= prim[XXX+2]*frac_flip;
      
      //Flip to carbon
      //prim[XXX+1] += prim[XXX]*frac_flip;
      //prim[XXX] -= prim[XXX]*frac_flip;
      
      for( q=0 ; q<NUM_I ; ++q ){ ebind2 += EBIND[q]*prim[q+XXX]; }
      prim[PPP] = etot + (ebind1-ebind2);

      prim[PPP] = get_pre_from_etot( prim , &temp );
      *T = temp;
   }
   
   /////////////////////////////////////////////////////
   /*
   int q;
   double etot = 0.0, ebind1 = 0.0, ebind2 = 0.0;
   etot = get_eint( prim , &temp );
   for( q=0 ; q<NUM_I ; ++q ){ ebind1 += EBIND[q]*prim[q+XXX]; }

   double flip = 0.30;//20% He4 to C12
   if( fabs(rr[l])<1e9 ){
      prim[XXX+1] += prim[XXX]*flip;
      prim[XXX] -= prim[XXX]*flip;
      
      for( q=0 ; q<NUM_I ; ++q ){ ebind2 += EBIND[q]*prim[q+XXX]; }
      prim[PPP] = etot + (ebind1-ebind2);

      prim[PPP] = get_pre_from_etot( prim , &temp );
      *T = temp;
   }
   */
   
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
