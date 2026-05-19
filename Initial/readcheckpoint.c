
#include "../paul.h"

void prim2cons( double * , double * , double , double , double * );
void cons2prim( double * , double * , double , double , double * );

static int NL = 0;
static double * rr  = NULL;
static double * dr  = NULL;
static double * rho = NULL;
static double * Pp  = NULL;
static double * vr  = NULL;
static double * ALP = NULL;
static double * XHE = NULL;
static double * XC  = NULL;
static double * XO  = NULL;
static double * XNE = NULL;
static double * XMG = NULL;
static double * XSI = NULL;
static double * XS  = NULL;
static double * XAR = NULL;
static double * XCA = NULL;
static double * XTI = NULL;
static double * XCR = NULL;
static double * XFE = NULL;
static double * XNI = NULL;
static double * Mas = NULL;
static double * ent = NULL;
static double * tem = NULL;
static double * cs  = NULL;
static double * Ein = NULL;

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
   int nL = countlines("Initial/checkpoint.dat");
   rr  = (double *) malloc( nL*sizeof(double) );
   dr  = (double *) malloc( nL*sizeof(double) );
   rho = (double *) malloc( nL*sizeof(double) );
   Pp  = (double *) malloc( nL*sizeof(double) );
   vr  = (double *) malloc( nL*sizeof(double) );
   ALP = (double *) malloc( nL*sizeof(double) );
   XHE = (double *) malloc( nL*sizeof(double) );
   XC  = (double *) malloc( nL*sizeof(double) );
   XO  = (double *) malloc( nL*sizeof(double) );
   XNE = (double *) malloc( nL*sizeof(double) );
   XMG = (double *) malloc( nL*sizeof(double) );
   XSI = (double *) malloc( nL*sizeof(double) );
   XS  = (double *) malloc( nL*sizeof(double) );
   XAR = (double *) malloc( nL*sizeof(double) );
   XCA = (double *) malloc( nL*sizeof(double) );
   XTI = (double *) malloc( nL*sizeof(double) );
   XCR = (double *) malloc( nL*sizeof(double) );
   XFE = (double *) malloc( nL*sizeof(double) );
   XNI = (double *) malloc( nL*sizeof(double) );
   Mas = (double *) malloc( nL*sizeof(double) );
   ent = (double *) malloc( nL*sizeof(double) );
   tem = (double *) malloc( nL*sizeof(double) );
   cs  = (double *) malloc( nL*sizeof(double) );
   Ein = (double *) malloc( nL*sizeof(double) );

   FILE * pFile = fopen("Initial/checkpoint.dat","r");
   int l;
   for( l=0 ; l<nL ; ++l ){
      fscanf(pFile,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le",&(rr[l]),&(dr[l]),&(rho[l]),&(Pp[l]),&(vr[l]),&(ALP[l]),&(XHE[l]),&(XC[l]),&(XO[l]),&(XNE[l]),&(XMG[l]),&(XSI[l]),&(XS[l]),&(XAR[l]),&(XCA[l]),&(XTI[l]),&(XCR[l]),&(XFE[l]),&(XNI[l]),&(Mas[l]),&(ent[l]),&(tem[l]),&(cs[l]),&(Ein[l])); 
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
   //if( l==0 ) ++l;
   l--;

   //double rp = rr[l];
   //double rm = rr[l-1];
  
   prim[RHO] = rho[l];
   prim[PPP] = Pp[l];
   prim[VRR] = vr[l];
   prim[AAA] = ALP[l];
   prim[XXX+0] = XHE[l];
   prim[XXX+1] = XC[l];
   prim[XXX+2] = XO[l];
   prim[XXX+3] = XNE[l];
   prim[XXX+4] = XMG[l];
   prim[XXX+5] = XSI[l];
   prim[XXX+6] = XS[l];
   prim[XXX+7] = XAR[l];
   prim[XXX+8] = XCA[l];
   prim[XXX+9] = XTI[l];
   prim[XXX+10] = XCR[l];
   prim[XXX+11] = XFE[l];
   prim[XXX+12] = XNI[l];

   double temp = tem[l];
   *T = temp;

}

void initial_grid( int ii , double * rr_i , double * dr_i ){
   *rr_i = rr[ii];
   *dr_i = dr[ii];
}

void freeTable( void ){
   free(rr);
   free(rho);
   free(Pp);
   free(vr);
   free(ALP);
   free(XHE);
   free(XC);
   free(XO);
   free(XNE);
   free(XMG);
   free(XSI);
   free(XS);
   free(XAR);
   free(XCA);
   free(XTI);
   free(XCR);
   free(XFE);
   free(XNI);
   free(Mas);
   free(ent);
   free(tem);
   free(cs);
   free(Ein);
}
