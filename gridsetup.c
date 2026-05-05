
#include "paul.h"

void setICparams( struct domain * );
void initial_grid( int , double * , double * );

int getN0( int drank , int dsize , int dnum ){
   int N0 = (dnum*drank)/dsize;
   return(N0);
}

void setupGrid( struct domain * theDomain ){

   setICparams( theDomain );

   int Ng = NUM_G;
   theDomain->Ng = Ng;
   int Num_R = theDomain->theParList.Num_R;
   int LogZoning = theDomain->theParList.LogZoning;

   double f = theDomain->theParList.LogFractn; //frac of zones in linear regime
   double n = theDomain->theParList.LogSmooth; //smoothness param
   double scal, Rlin, Rlog;

   int rank = theDomain->rank;
   int size = theDomain->size;

   double Rmin = theDomain->theParList.rmin;
   double Rmax = theDomain->theParList.rmax;

   int N0r = getN0( rank   , size , Num_R );
   int N1r = getN0( rank+1 , size , Num_R );
   if( rank != 0 ) N0r -= Ng;
   if( rank != size-1 ) N1r += Ng;
   int Nr = N1r-N0r;

   theDomain->Nr = Nr;
   theDomain->theCells = (struct cell *) malloc( Nr*sizeof(struct cell));
   printf("Rank = %d, Nr = %d\n",theDomain->rank,Nr);

   int i;

   double dx = 1./(double)Num_R;
   double x0 = (double)N0r/(double)Num_R;
   double R0 = theDomain->theParList.LogRadius;

   double rr, dr;
   if( LogZoning == 4 ){
      //Read position in from checkpoint
      for( i=0 ; i<Nr ; ++i ){
         initial_grid( i+N0r , &rr , &dr );
         theDomain->theCells[i].riph = rr + 0.5*dr;
         theDomain->theCells[i].dr   = dr;
      }
   }else{
      //Log/Linear position distribution
      for( i=0 ; i<Nr ; ++i ){
         double xm = x0 + ((double)i   )*dx;
         double xp = x0 + ((double)i+1.)*dx;
         double rp,rm;
         if( LogZoning == 0 ){
            rp = Rmin + xp*(Rmax-Rmin);
            rm = Rmin + xm*(Rmax-Rmin);
            if( (rank == 0) & (i==0) ){ rm = 0.0; }
         }else if( LogZoning == 1 ){
            rp = Rmin*pow(Rmax/Rmin,xp);
            rm = Rmin*pow(Rmax/Rmin,xm);
         }else if( LogZoning == 2 ){
            //rp = R0*pow(Rmax/R0,xp) + Rmin-R0 + (R0-Rmin)*xp;
            //rm = R0*pow(Rmax/R0,xm) + Rmin-R0 + (R0-Rmin)*xm;
            rp = (Rmax-Rmin)*(pow(Rmax/R0,xp)-1.0)/(Rmax/R0-1.0) + Rmin;
            rm = (Rmax-Rmin)*(pow(Rmax/R0,xm)-1.0)/(Rmax/R0-1.0) + Rmin;
         }else{ //LogZoning == 3
            scal = pow(pow((R0-Rmin)/f,n)+pow(Rmax-R0*pow(Rmax/R0,-f/(1.-f)),n),1./n);
            //Outer radius// 
            Rlin = xp*(R0-Rmin)/f;
            Rlog = R0*pow(Rmax/R0,(xp-f)/(1.-f)) - R0*pow(Rmax/R0,-f/(1.-f));
            rp = (Rmax-Rmin)*pow(pow(Rlin,n)+pow(Rlog,n),1./n)/scal + Rmin;
            //Inner radius// 
            Rlin = xm*(R0-Rmin)/f;
            Rlog = R0*pow(Rmax/R0,(xm-f)/(1.-f)) - R0*pow(Rmax/R0,-f/(1.-f));
            rm = (Rmax-Rmin)*pow(pow(Rlin,n)+pow(Rlog,n),1./n)/scal + Rmin;
         }
         theDomain->theCells[i].riph = rp;
         theDomain->theCells[i].dr   = rp - rm;
      }
   }

}


