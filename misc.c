#include "paul.h"
#include <string.h>

double get_dA( double );
double get_dV( double , double );
double get_g( struct cell * );
double get_GMr( struct cell * );
double get_moment_arm( double , double );

void cons2prim( double * , double * , double , double , double * );

double mindt( double * , double , double , double , double , double * );

double getmindt( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Ng = theDomain->Ng;
   int gE = theDomain->theParList.grav_e_mode;
   int gFl = theDomain->theParList.grav_flag;

   double dt = 1e100;
   int imin = Ng;
   int imax = Nr-Ng;
   if( theDomain->rank == 0) imin = 0;
   if( theDomain->rank == theDomain->size-1 ) imax = Nr;

   int i;
   for( i=imin ; i<imax ; ++i ){
      int im = i-1;
      if( i==0 ) im = 0;
      struct cell * c = theCells+i;
      double dr = c->dr;
      double r = c->riph-.5*dr;
      double wm = theCells[im].wiph;
      double wp = c->wiph;
      double w = .5*(wm+wp);
      double g = 0.0;
      if( gFl && gE == 1 ) g = get_g( c );
      double dt_temp = mindt( c->prim , w , r , g , dr , &c->T );
      if( dt > dt_temp ) dt = dt_temp;
   }
   dt *= theDomain->theParList.CFL; 
   MPI_Allreduce( MPI_IN_PLACE , &dt , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD );

   return( dt );
}

double get_vr( double * );

void set_wcell( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int mesh_motion = theDomain->theParList.Mesh_Motion;
   int Nr = theDomain->Nr;
   int bufferzone = 1;
   double Rmax = theCells[Nr-1].riph;
   MPI_Allreduce( MPI_IN_PLACE , &Rmax , 1 , MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD );
   double Rbuf = 0.8*Rmax;

   int i;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * cL = theCells+i;  
      double w = 0.0;
      if( mesh_motion && i<Nr-1 ){
         struct cell * cR = theCells+i+1;
         double wL = get_vr( cL->prim );
         double wR = get_vr( cR->prim );
         w = .5*(wL + wR); 
         //if( i==0 && theDomain->rank==0 ) w = 0.5*wR*(cL->riph)/(cR->riph-.5*cR->dr);//*(cR->riph - .5*cR->dr)/(cL->riph);//0.0;//2./3.*wR;
         if( bufferzone && cL->riph > Rbuf ) w *= (Rmax-cL->riph)/(Rmax-Rbuf);//w = 0.0;
         //if( i==0 && theDomain->rank==0 ) w = 0.0;
      }
      cL->wiph = w;
   }
}

void adjust_RK_cons( struct domain * theDomain , double RK ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   int i,q;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      for( q=0 ; q<NUM_Q ; ++q ){
         c->cons[q] = (1.-RK)*c->cons[q] + RK*c->RKcons[q];
      }
   }
}

void move_cells( struct domain * theDomain , double RK , double dt){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int i;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      c->riph += c->wiph*dt;
   }
}

void calc_dr( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   int i;
   for( i=1 ; i<Nr ; ++i ){
      int im = i-1;
      double rm = theCells[im].riph;
      double rp = theCells[i ].riph;
      double dr = rp-rm;
      theCells[i].dr = dr;
   }
   if( theDomain->rank==0 ) theCells[0].dr = theCells[0].riph;
}

void calculate_mass( struct domain * );

void calc_prim( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int gE = theDomain->theParList.grav_e_mode;
   if( gE ) calculate_mass( theDomain );

   int Nr = theDomain->Nr;
   int Ng = theDomain->Ng;
   int imin = Ng;
   int imax = Nr-Ng;
   if( theDomain->rank==0 ){ imin=0; }
   if( theDomain->rank==theDomain->size-1 ){ imax=Nr; }
   int i;
   for( i=imin ; i<imax ; ++i ){
      struct cell * c = theCells+i;
      double rp = c->riph;
      double rm = rp-c->dr;
      double dV = get_dV( rp , rm );
//      double g = 0.0;
//      if( gE == 1 ) g = get_g( c );
      double GMr = 0.0;
      if( gE == 3 ) GMr = get_GMr( c );
      cons2prim( c->cons , c->prim , GMr , dV , &c->T );
   }
}

void plm( struct domain * );
void riemann( struct cell * , struct cell * , double , double );
void calculate_fgrav( struct domain * );

void radial_flux( struct domain * theDomain , double dt ){

   int gE = theDomain->theParList.grav_e_mode;
   if( gE ) calculate_mass(theDomain);
   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   plm( theDomain );

   int i;
   int imin = 1;
   if( theDomain->rank==0 ){ imin=0; }
   for( i=imin ; i<Nr-1 ; ++i ){
      struct cell * cL = theCells+i;
      struct cell * cR = theCells+i+1;
      double r = cL->riph;
      double dA = get_dA(r); 
      riemann( cL , cR , r , dA*dt );
   }
}

void source( double * , double * , double , double , double );
void source_nozz( double * , double * , double , double , double , double );
void source_alpha( double * , double * , double * , double , double , double * );
void source_grow( double * , double * , double * , double , double );
void gravity_addsrc( struct domain * , double );

void add_source( struct domain * theDomain , double dt ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   double t = theDomain->t;
   double grad[NUM_Q];

   int i,q;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      double rp = c->riph;
      double rm = rp-c->dr;
      double r = get_moment_arm(rp,rm);
      double dV = get_dV(rp,rm);
      source( c->prim , c->cons , rp , rm , dV*dt );
      source_nozz( c->prim , c->cons , rp , rm , t , dV*dt );
      int inside = i>0 && i<Nr-1;
      for( q=0 ; q<NUM_Q ; ++q ){
         if( inside ){
            struct cell * cp = theCells+i+1;
            struct cell * cm = theCells+i-1;
            double dR = .5*cp->dr + c->dr + .5*cm->dr;
            grad[q] = (cp->prim[q]-cm->prim[q])/dR;
/*
            double gR = 2.*(cp->prim[q]-c->prim[q])/(cp->dr+c->dr);
            double gL = 2.*(c->prim[q]-cm->prim[q])/(c->dr+cm->dr);
            double g = gL;
            if( gL*gR < 0. ) g = 0.;
            if( fabs(gR) < fabs(g) ) g = gR;
            grad[q] = g;
*/
         }else{
            grad[q] = 0.0;
         }
      }
      source_alpha( c->prim , c->cons , grad , r , dV*dt , &c->T );
      if( 0 ) source_grow( c->prim , c->cons , grad , r , dV*dt );
   }   

   int gravity_flag = theDomain->theParList.grav_flag;
   if( gravity_flag ) gravity_addsrc( theDomain , dt );
}


void longandshort( struct domain * theDomain , double * L , double * S , int * iL , int * iS , int * rL , int * rS ){ 

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   double rmax = theCells[Nr-1].riph;
   double rmin = theCells[0].riph;
   double R0 = theDomain->theParList.LogRadius;
   double f = theDomain->theParList.LogFractn;
   double Rmin0 = theDomain->theParList.rmin;
   MPI_Allreduce( MPI_IN_PLACE , &rmax , 1 , MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &rmin , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD );
   int Nr0 = theDomain->theParList.Num_R;
   double dr0 = rmax/(double)Nr0;
   double dx0 = log(rmax/rmin)/Nr0;
   int logscale = theDomain->theParList.LogZoning;

   double Long  = 0.0; 
   double Short = 0.0; 
   int iLong  = -1;
   int iShort = -1;

   int rank = theDomain->rank;
   int size = theDomain->size;
   int Ng   = theDomain->Ng;

   int imin = 0;
   if( logscale==1 ) imin=1;
   int imax = Nr;
   if( rank!=0 )      imin = Ng;
   if( rank!=size-1 ) imax = Nr-Ng;

   int i;
   for( i=imin ; i<imax ; ++i ){
      struct cell * c = theCells+i;
      double dy = c->dr;
      double dxS = dr0;
      double dxL = dr0;
      if( logscale==1 ) {
         dxS = c->riph*dx0; 
         dxL = c->riph*dx0;
      }
      if( logscale==2 ) {
         dxS = (1./(double)Nr0)*( c->riph-rmin + rmax/(rmax/R0-1.) )*log(rmax/R0);
         dxL = (1./(double)Nr0)*( c->riph-rmin + rmax/(rmax/R0-1.) )*log(rmax/R0);
      }
      if( logscale==3 ) {
         dy = c->dr/c->riph;
         dxS = fmin((R0-Rmin0)/(R0*f*(double)Nr0),(R0-Rmin0)/(c->riph*f*(double)Nr0));
         dxL = fmax((R0-Rmin0)/(R0*f*(double)Nr0),(R0-Rmin0)/(c->riph*f*(double)Nr0));
      }
      double l = dy/dxL;
      double s = dxS/dy;
      if( Long  < l ){ Long  = l; iLong  = i; } 
      if( Short < s ){ Short = s; iShort = i; } 
   }

   struct { double value ; int index ; } maxminbuf;
   maxminbuf.value = Short;
   maxminbuf.index = rank;
   MPI_Allreduce( MPI_IN_PLACE , &maxminbuf , 1 , MPI_DOUBLE_INT , MPI_MAXLOC , MPI_COMM_WORLD );
   *S = maxminbuf.value;
   *rS = maxminbuf.index;

   maxminbuf.value = Long;
   maxminbuf.index = rank;
   MPI_Allreduce( MPI_IN_PLACE , &maxminbuf , 1 , MPI_DOUBLE_INT , MPI_MAXLOC , MPI_COMM_WORLD );
   *L = Long;
   *rL = maxminbuf.index;

   *iS = iShort;
   *iL = iLong;

}

void shortoverboundary( struct domain * theDomain , int iS , int rS , int * bdry ){
   // 0 = short as normal ; 1 = short across left bdry ; 2 = short across right bdry
   int bdry_flag = 0;

   int rank = theDomain->rank;
   int size = theDomain->size;
   struct cell * theCells = theDomain->theCells;
   int Ng = theDomain->Ng;
   int Nr = theDomain->Nr;

   if ( rank == rS ){
      int iSp = iS+1;
      int iSm = iS-1;
      if( rank == size-1 && iS == Nr-1 ) iSp=iS;
      if( rank == 0 && iS == 0 ) iSm=0;

      int imin = Ng;
      int imax = Nr-Ng;
      if( rank==0 ) imin = 0;
      if( rank==size-1 ) imax = Nr;

      double drL = theCells[iSm].dr;
      double drR = theCells[iSp].dr;
      if( drL<drR && iS==imin ) bdry_flag = 1;
      if( drR<drL && iS==imax-1 ) bdry_flag = 2; //imax-1 if you want it to happen
   }

   MPI_Allreduce( MPI_IN_PLACE , &bdry_flag , 1 , MPI_INT , MPI_MAX , MPI_COMM_WORLD );
   *bdry = bdry_flag;
}

void AMR( struct domain * theDomain ){

   double L,S;
   int iL=0;
   int iS=0;
   int rL=0;
   int rS=0;
   longandshort( theDomain , &L , &S , &iL , &iS , &rL , &rS );
   int bdry;
   shortoverboundary( theDomain , iS , rS , &bdry );
   int rank = theDomain->rank;
   int size = theDomain->size;

   double MaxShort = theDomain->theParList.MaxShort;
   double MaxLong  = theDomain->theParList.MaxLong;

   struct cell * theCells = theDomain->theCells;
   int Ng = theDomain->Ng;
   int Nr = theDomain->Nr;

   int gE = theDomain->theParList.grav_e_mode;

   if( S>MaxShort && MaxShort != 0 ){
      ////////////////////// Normal Shorting ///////////////////////////
      if( bdry == 0 && rank == rS ){
         printf("KILL! Rank %d; Short = %e #%d of %d\n",rank,S,iS,Nr);

         int iSp = iS+1;
         int iSm = iS-1;
         if( rank == size-1 && iS == Nr-1 ) iSp=iS;
         if( rank == 0 && iS == 0 ) iSm=0;
         //Possibly shift iS backwards by 1 
         double drL = theCells[iSm].dr;
         double drR = theCells[iSp].dr;
         if( drL<drR ){
            --iS;
            --iSm;
            --iSp;
         }
         struct cell * c  = theCells+iS;
         struct cell * cp = theCells+iSp;

         //Remove Zone at iS+1
         c->dr   += cp->dr;
         c->riph  = cp->riph;
         c->dm   += cp->dm;
         c->miph  = cp->miph;
         int q;
         for( q=0 ; q<NUM_Q ; ++q ){
            c->cons[q]   += cp->cons[q];
            c->RKcons[q] += cp->RKcons[q];
         }
         double rp = c->riph;
         double rm = rp - c->dr;
         double dV = get_dV( rp , rm );
   //      double g = 0.0;
   //      if( gE == 1 ) g = get_g( c );
         double GMr = 0.0;
         if( gE == 3 ) GMr = get_GMr( c );
         cons2prim( c->cons , c->prim , GMr , dV , &c->T );

         //Shift Memory
         int blocksize = Nr-iSp-1;
         memmove( theCells+iSp , theCells+iSp+1 , blocksize*sizeof(struct cell) );
         theDomain->Nr -= 1;
         Nr = theDomain->Nr;
         theDomain->theCells = (struct cell *) realloc( theCells , Nr*sizeof(struct cell) );
         theCells = theDomain->theCells;
         if( iS < iL ) iL--;
      }

      ////////////////////// Shorting over Left ///////////////////////////
      if( bdry == 1 && rank == rS ){
         printf("KILL! LEFT Rank %d; Short = %e #%d of %d\n",rank,S,iS,Nr);

         int blocksize = Nr-Ng-1;
         memmove( theCells+Ng , theCells+Ng+1 , blocksize*sizeof(struct cell) );
         theDomain->Nr -= 1;
         Nr = theDomain->Nr;
         theDomain->theCells = (struct cell *) realloc( theCells , Nr*sizeof(struct cell) );
         theCells = theDomain->theCells;
         iL--;
      }
      if( bdry == 1 && rank == rS-1 ){
         printf("KILL! LEFT Rank %d; Merge #%d and #%d\n",rank,Nr-Ng-1,Nr-Ng);

         struct cell * c  = theCells+Nr-Ng-1;
         struct cell * cp = theCells+Nr-Ng;

         c->dr   += cp->dr;
         c->riph  = cp->riph;
         c->dm   += cp->dm;
         c->miph  = cp->miph;
         int q;
         for( q=0 ; q<NUM_Q ; ++q ){
            c->cons[q]   += cp->cons[q];
            c->RKcons[q] += cp->RKcons[q];
         }
         double rp = c->riph;
         double rm = rp - c->dr;
         double dV = get_dV( rp , rm );
   //      double g = 0.0;
   //      if( gE == 1 ) g = get_g( c );
         double GMr = 0.0;
         if( gE == 3 ) GMr = get_GMr( c );
         cons2prim( c->cons , c->prim , GMr , dV , &c->T );
      }
      ////////////////////// Shorting over Right ///////////////////////////
      if( bdry == 2 && rank == rS+1 ){
         printf("KILL! RIGHT Rank %d; Remove Zone #%d of #%d\n",rank,Ng,Nr);

         int blocksize = Nr-Ng-1;
         memmove( theCells+Ng , theCells+Ng+1 , blocksize*sizeof(struct cell) );
         theDomain->Nr -= 1;
         Nr = theDomain->Nr;
         theDomain->theCells = (struct cell *) realloc( theCells , Nr*sizeof(struct cell) );
         theCells = theDomain->theCells;
         iL--;
      }
      if( bdry == 2 && rank == rS ){
         printf("KILL! RIGHT Rank %d; Short = %e #%d of %d\n",rank,S,iS,Nr);
         
         struct cell * c  = theCells+Nr-Ng-1;
         struct cell * cp = theCells+Nr-Ng;

         c->dr   += cp->dr;
         c->riph  = cp->riph;
         c->dm   += cp->dm;
         c->miph  = cp->miph;
         int q;
         for( q=0 ; q<NUM_Q ; ++q ){
            c->cons[q]   += cp->cons[q];
            c->RKcons[q] += cp->RKcons[q];
         }
         double rp = c->riph;
         double rm = rp - c->dr;
         double dV = get_dV( rp , rm );
   //      double g = 0.0;
   //      if( gE == 1 ) g = get_g( c );
         double GMr = 0.0;
         if( gE == 3 ) GMr = get_GMr( c );
         cons2prim( c->cons , c->prim , GMr , dV , &c->T );
      }
   }

   if( L>MaxLong && MaxLong != 0 ){
      if( rank==rL ){
         printf("FORGE! Rank %d; Long  = %e #%d of %d\n",rank,L,iL,Nr);

         theDomain->Nr += 1;
         Nr = theDomain->Nr;
         theDomain->theCells = (struct cell *) realloc( theCells , Nr*sizeof(struct cell) );
         theCells = theDomain->theCells;
         int blocksize = Nr-iL-1;
         memmove( theCells+iL+1 , theCells+iL , blocksize*sizeof(struct cell) );

         struct cell * c  = theCells+iL;
         struct cell * cp = theCells+iL+1;

         double rp = c->riph;
         double rm = rp - c->dr;
         double r0 = pow( .5*(rp*rp*rp+rm*rm*rm) , 1./3. );
         double dm = .5*c->dm;

         c->riph  = r0;
         c->dr    = r0-rm;
   //      cp->riph = rp;
         cp->dr   = rp-r0;

         c->dm    = dm;
         cp->dm   = dm;
   //      cp->miph = c->miph;
         c->miph -= dm;

         int q;
         for( q=0 ; q<NUM_Q ; ++q ){
            c->cons[q]    *= .5;
            c->RKcons[q]  *= .5;
            cp->cons[q]   *= .5;
            cp->RKcons[q] *= .5;
         }

         double dV = get_dV( r0 , rm );
   //      double g = 0.0;
         double GMr = 0.0;
   //      if( gE == 1 ) g = get_g( c );
         if( gE == 3 ) GMr = get_GMr( c );
         cons2prim( c->cons , c->prim , GMr , dV , &c->T );
         dV = get_dV( rp , r0 );
   //      if( gE == 1 ) g = get_g( cp );
         if( gE == 3 ) GMr = get_GMr( cp );
         cons2prim( cp->cons , cp->prim , GMr , dV , &cp->T );
      }
   }
}


