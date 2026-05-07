
#include "../../paul.h"

//functions from euler-'EOS'.c
void init_eos( int );
void prim2cons( double * , double * , double , double , double * );
void cons2prim( double * , double * , double , double , double * );

//functions from geometry.c
double get_dV( double , double );

static double * rr  = NULL;
static double * dr  = NULL;
static double * rho = NULL;
static double * Pp  = NULL;
static double * vr  = NULL;
static double * XXA = NULL;
static double * XX0 = NULL;
static double * XX1 = NULL;
static double * XX2 = NULL;
static double * XX3 = NULL;
static double * XX4 = NULL;
static double * XX5 = NULL;
static double * XX6 = NULL;
static double * XX7 = NULL;
static double * XX8 = NULL;
static double * XX9 = NULL;
static double * X10 = NULL;
static double * X11 = NULL;
static double * X12 = NULL;
static double * miph= NULL;
static double * stot= NULL;
static double * Temp= NULL;
static double * cs  = NULL;
static double * etot= NULL;

static double * logr= NULL;
static double * logv= NULL;
static double * logs= NULL;

int main( int argc , char * argv[] ){
   init_eos(2);

   /////////////////////////////////////////////////////////////////////////
   //                              INPUTS                                 //
   char filename[] = "output.dat";
   
   //           Crop parameters : homologous expansion fit                //
   double crop_r_lo = 1e6;
   double crop_r_hi = 5.9e10;

   int nL_sedona = 137;
   char outfilename[] = "initial.mod";

   /////////////////////////////////////////////////////////////////////////
   
   //count number of hydro zones//
   FILE * pFile = fopen(filename, "r");
   int nL=0;
   char c;
   while ((c = fgetc(pFile)) != EOF){
      if (c == '\n') ++nL;
   }
   fclose(pFile);
   nL-=1; //dont count the header

   printf("Begin reading input file\n\n");

   FILE * fp;
   char header_buffer[500];

   fp = fopen(filename, "r");

   if( fp == NULL ){
      perror("Error opening file");
      return 1;
   }

   // Read and discard the entire header line using fgets
   // Check for NULL return to handle errors or end of file
   if( fgets(header_buffer, sizeof(header_buffer), fp) == NULL ){
      perror("Error reading header line");
      fclose(fp);
      return 1;
   }

   int error;
   rr  = (double *) malloc( nL*sizeof(double) );
   dr  = (double *) malloc( nL*sizeof(double) );
   rho = (double *) malloc( nL*sizeof(double) );
   Pp  = (double *) malloc( nL*sizeof(double) );
   vr  = (double *) malloc( nL*sizeof(double) );
   XXA = (double *) malloc( nL*sizeof(double) );
   XX0 = (double *) malloc( nL*sizeof(double) );
   XX1 = (double *) malloc( nL*sizeof(double) );
   XX2 = (double *) malloc( nL*sizeof(double) );
   XX3 = (double *) malloc( nL*sizeof(double) );
   XX4 = (double *) malloc( nL*sizeof(double) );
   XX5 = (double *) malloc( nL*sizeof(double) );
   XX6 = (double *) malloc( nL*sizeof(double) );
   XX7 = (double *) malloc( nL*sizeof(double) );
   XX8 = (double *) malloc( nL*sizeof(double) );
   XX9 = (double *) malloc( nL*sizeof(double) );
   X10 = (double *) malloc( nL*sizeof(double) );
   X11 = (double *) malloc( nL*sizeof(double) );
   X12 = (double *) malloc( nL*sizeof(double) );
   miph= (double *) malloc( nL*sizeof(double) );
   stot= (double *) malloc( nL*sizeof(double) );
   Temp= (double *) malloc( nL*sizeof(double) );
   cs  = (double *) malloc( nL*sizeof(double) );
   etot= (double *) malloc( nL*sizeof(double) );

   for( int l=0 ; l<nL ; l++){
   // Read in the data //
      error = fscanf(fp,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",\
         &(rr[l]),&(dr[l]),\
         &(rho[l]),&(Pp[l]),&(vr[l]),&(XXA[l]),\
         &(XX0[l]),&(XX1[l]),&(XX2[l]),&(XX3[l]),&(XX4[l]),\
         &(XX5[l]),&(XX6[l]),&(XX7[l]),&(XX8[l]),&(XX9[l]),\
         &(X10[l]),&(X11[l]),&(X12[l]),\
         &(miph[l]),&(stot[l]),&(Temp[l]),&(cs[l]),&(etot[l]));
   }

   fclose(fp);

   /////////////////////////////////////////////////////////////////////////
   //              Linear Regression in log log space                     //
   printf("Begin Linear Regression in log log space\n");
   int i_lo=-1, i_hi=-1;

   for( int l=0 ; l<nL ; l++){
      if( rr[l] > crop_r_lo && i_lo < 0 ){ i_lo = l; }
      if( rr[l] > crop_r_hi && i_hi < 0 ){ i_hi = l; }
   }

   int nL_trimmed = i_hi-i_lo+1;

   //copy radius and velocity arrays
   logr= (double *) malloc( nL*sizeof(double) );
   logv= (double *) malloc( nL*sizeof(double) );
   logs= (double *) malloc( nL*sizeof(double) );
   memcpy(logr, rr, nL*sizeof(double));
   memcpy(logv, vr, nL*sizeof(double));
   memcpy(logs, rr, nL*sizeof(double));
   //shift the arrays down
   memmove(logr, &logr[i_lo], nL_trimmed*sizeof(double));
   memmove(logv, &logv[i_lo], nL_trimmed*sizeof(double));
   memmove(logs, &logs[i_lo], nL_trimmed*sizeof(double));
   //trim the arrays down
   double * temp; 
   temp = realloc(logr, nL_trimmed*sizeof(double));
   logr = temp;
   temp = realloc(logv, nL_trimmed*sizeof(double));
   logv = temp;
   temp = realloc(logs, nL_trimmed*sizeof(double));
   logs = temp;

   for( int l=0 ; l<nL_trimmed ; l++){
      logr[l] = log10(logr[l]);
      logv[l] = log10(logv[l]);
      //logs[l] = log10(logs[l]);
      logs[l] = 1.0;
   }
   
   // Numerical Recipes in C -- Chap 15.2 //
   //x = logr , y = logr , sig = 1/weight^2 = logs , log(vr) = P1*log(rr) + P0
   double wt, t, sxoss, sx=0, sy=0, st2=0, ss=0, P1=0, P0;

   for( int i=0 ; i<nL_trimmed ; i++ ){ //...with weights
      wt=1.0/sqrt(logs[i]);
      ss += wt;
      sx += logr[i]*wt;
      sy += logv[i]*wt;
   }
   sxoss=sx/ss;
   for( int i=0 ; i<nL_trimmed ; i++ ){
      t=(logr[i]-sxoss)/logs[i];
      st2 += t*t;
      P1 += t*logv[i]/logs[i];
   }
   P1 /= st2;
   P0 = (sy-sx*P1)/ss;
   
   double t_exp = pow(10.,-P0/P1);
   printf("Homologous Expansion Fit: P1, P0 = %e, %e\n",P1,P0);
   printf("Time since explosion: t_exp = %e\n\n",t_exp);

   /////////////////////////////////////////////////////////////////////////
   double Mtot = miph[nL-1];
   double mass_enc;
   double bucket_mass = Mtot/(double)nL_sedona;

   FILE * fpp = fopen( outfilename , "w" );
   fprintf(fpp,"1D_sphere standard\n");
   fprintf(fpp,"%d 0.0 %e 13\n",nL_sedona-1,t_exp);
   fprintf(fpp,"2.4 6.12 8.16 10.20 12.24 14.28 16.32 18.36 20.40 22.44 24.48 26.52 28.56\n");
   fclose(fpp);

   double rp, rm, dV, GMr, T, rp_bin, rm_bin=0.0;
   double M,r2_3,r4_5,mt;
   double prim[NUM_Q],cons[NUM_Q],prim_bin[NUM_Q],cons_bin[NUM_Q];
   for( int i=0 ; i<NUM_Q ; ++i ){ cons_bin[i] = 0; }
   int left = 0;
   int ilog = 0;
   for( int l=0 ; l<nL ; l++){
      FILE * fpp = fopen( outfilename , "a" );

      //get volume//
      rp = rr[l] + 0.5*dr[l];
      rm = rp - dr[l];
      dV = get_dV( rp , rm );

      //get GMr//
      M = miph[l] - rho[l]*dV;
      r2_3 = (rp*rp + rm*rm + rp*rm)/3.;
      r4_5 = (pow(rp,4.) + pow(rp,3.)*rm + rp*rp*rm*rm + rp*pow(rm,3.) + pow(rm,4.) )/5.;
      mt = M - 4./3.*M_PI*rm*rm*rm*rho[l];
      GMr = grav_G/r2_3*( mt*rr[l] + 4./3.*M_PI*rho[l]*r4_5 );

      //prim2cons
      prim[RHO]      = rho[l];
      prim[PPP]      = Pp[l];
      prim[VRR]      = vr[l];
      prim[AAA]      = XXA[l];
      prim[XXX]      = XX0[l];
      prim[XXX+1]    = XX1[l];
      prim[XXX+2]    = XX2[l];
      prim[XXX+3]    = XX3[l];
      prim[XXX+4]    = XX4[l];
      prim[XXX+5]    = XX5[l];
      prim[XXX+6]    = XX6[l];
      prim[XXX+7]    = XX7[l];
      prim[XXX+8]    = XX8[l];
      prim[XXX+9]    = XX9[l];
      prim[XXX+10]   = X10[l];
      prim[XXX+11]   = X11[l];
      prim[XXX+12]   = X12[l];
      T = Temp[l];

      prim2cons(prim, cons, GMr, dV, &T);

      //fill the bucket//
      for( int i=0 ; i<NUM_Q ; ++i ){ cons_bin[i] += cons[i]; }

      //check if we reached mass bucket//
      if( cons_bin[DDD] >= bucket_mass ){  //skip final bucket
      //if( cons_bin[DDD] >= bucket_mass || l == nL-1 ){  //write out final bucket
         printf("writing line %d\n",l);

         //bucket_mass = Mtot/log(10)*((10.-pow(10.,1./(double)nL_sedona))/((double)nL_sedona-1.));
         //bucket_mass /= ((double)ilog/((double)nL_sedona-1.))*(10.-pow(10.,1./(double)nL_sedona)) + pow(10.,1./(double)nL_sedona);
         printf("Bucket size %e\n",bucket_mass);
         //bucket_mass = Mtot/(double)nL_sedona;
         
         bucket_mass = -log10((double)ilog/((double)nL_sedona)*(10.-pow(10.,1./(double)nL_sedona))+pow(10.,1./(double)nL_sedona));
         ilog++;
         bucket_mass += log10((double)ilog/((double)nL_sedona)*(10.-pow(10.,1./(double)nL_sedona))+pow(10.,1./(double)nL_sedona));
         bucket_mass *= Mtot;

         //get volume//
         rp_bin = rr[l] + 0.5*dr[l];
         rm_bin = rr[left] - 0.5*dr[left];
         dV = get_dV( rp_bin , rm_bin );

         //get GMr//
         M = miph[left-1];
         r2_3 = (rp*rp + rm*rm + rp*rm)/3.;
         r4_5 = (pow(rp,4.) + pow(rp,3.)*rm + rp*rp*rm*rm + rp*pow(rm,3.) + pow(rm,4.) )/5.;
         mt = M - 4./3.*M_PI*rm*rm*rm*cons_bin[DDD]/dV;
         GMr = grav_G/r2_3*( mt*0.5*(rp_bin+rm_bin) + 4./3.*M_PI*cons_bin[DDD]/dV*r4_5 );

         cons2prim(cons_bin, prim_bin, GMr, dV, &T);

         //////////////////////////////////////////////////////////
         //Hydro Comparison
         //fprintf(fpp,"%.14e %.14e %.14e %.14e ",0.5*(rp_bin+rm_bin),prim_bin[VRR],prim_bin[RHO],T);
         
         //Velocity is local fluid velocity
         fprintf(fpp,"%.14e %.14e %.14e %.14e ",rp_bin,vr[l],prim_bin[RHO],T);
         
         //Velocity is homologous -- no temp recalc
         //fprintf(fpp,"%.14e %.14e %.14e %.14e ",rp_bin,rp_bin/t_exp,prim_bin[RHO],T);
         //////////////////////////////////////////////////////////
         
         for( int i=0 ; i<NUM_I ; ++i ){ fprintf(fpp,"%.14e ",prim_bin[XXX+i]); }
         fprintf(fpp,"\n");

         //empty the bucket//
         printf("Empty %e grams\n",cons_bin[DDD]);
         for( int i=0 ; i<NUM_Q ; ++i ){ cons_bin[i] = 0; }
         left = l+1;
      }

      fclose(fpp);
   }

   /////////////////////////////////////////////////////////////////////////   

   free(rr);
   free(dr);
   free(rho);
   free(Pp);
   free(vr);
   free(XXA);
   free(XX0);
   free(XX1);
   free(XX2);
   free(XX3);
   free(XX4);
   free(XX5);
   free(XX6);
   free(XX7);
   free(XX8);
   free(XX9);
   free(X10);
   free(X11);
   free(X12);
   free(miph);
   free(stot);
   free(Temp);
   free(cs);
   free(etot);

   free(logr);
   free(logv);
   free(logs);
   return(0);

}

