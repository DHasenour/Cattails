
#include "paul.h"

static int grav_E_mode = 0;
static double grav_G = 0.0;

void setGravityParams( struct domain * theDomain ){
   grav_E_mode = theDomain->theParList.grav_e_mode;
   grav_G = theDomain->theParList.grav_G;
}

double get_dV( double , double );

double get_g( struct cell * c ){

   double G = grav_G;
   double rp = c->riph;
   double rm = c->riph - c->dr;


   double Mm = c->miph - c->dm;
   double rho = c->dm/(4./3.*M_PI*(rp*rp*rp-rm*rm*rm));

   double r2_3 = (rp*rp + rm*rm + rp*rm)/3.;
   double r3_4 = (rp*rp*rp + rp*rp*rm + rp*rm*rm + rm*rm*rm)/4.;

   double g = -G/r2_3*( Mm + 4./3.*M_PI*rho*(r3_4 - rm*rm*rm) );

   return( g );
 
}

double get_GMr( struct cell * c ){
   double G = grav_G;
   double rp = c->riph;
   double rm = rp - c->dr;
   double rc = .5*(rp+rm);
   double M = c->miph - c->dm;
   double dV = 4./3.*M_PI*(rp*rp*rp-rm*rm*rm);
   double rho = c->dm/dV;
//   double r = (rp*rp+rm*rm+rp*rm)/3./rc;
   double r2_3 = (rp*rp + rm*rm + rp*rm)/3.;
   double r4_5 = (pow(rp,4.) + pow(rp,3.)*rm + rp*rp*rm*rm + rp*pow(rm,3.) + pow(rm,4.) )/5.;

   double mt = M - 4./3.*M_PI*rm*rm*rm*rho;

   double eps = G/r2_3*( mt*rc + 4./3.*M_PI*rho*r4_5 );

   return( eps );
}

void aggregate_mass( struct domain * theDomain ){

   int rank = theDomain->rank;
   int size = theDomain->size;
   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Ng = theDomain->Ng;

   double M = theDomain->point_mass;
   int i;
   int imin=Ng;
   if( rank==0 ) imin = 0;

   double Mtemp = M;
   for( i=imin-1 ; i>=0 ; --i ){
      theCells[i].miph = Mtemp;
      //Mtemp -= theCells[i+1].dm; //BUG ALERT//
      Mtemp -= theCells[i].dm;
   }

   for( i=imin ; i<Nr ; ++i ){
      //if( i>=imin ) M += theCells[i].dm;
      M += theCells[i].dm;
      theCells[i].miph = M;
   }

   int imax=Nr-Ng;
   if( rank==size-1 ) imax = Nr;
   double Mtot = theCells[imax-1].miph;
   double Mtot_inner = 0.0;
   double Mrecv = 0.0;

   int nrk;
   for( nrk=0 ; nrk < size ; ++nrk ){
      if( rank==nrk ){
         Mtot_inner = Mrecv;
         Mtot += Mtot_inner;
         if( nrk<size-1 ){
            MPI_Send( &Mtot , 1 , MPI_DOUBLE , nrk+1 , 666 , MPI_COMM_WORLD );
         }
      }
      if( rank==nrk+1 ){
         MPI_Status status;
         MPI_Recv( &Mrecv , 1 , MPI_DOUBLE , nrk , 666 , MPI_COMM_WORLD , &status );
      }
   }

   for( i=0 ; i<Nr ; ++i ){
      theCells[i].miph += Mtot_inner;
   }

}

void calculate_mass( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   int i;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      //double rp,rm;
      //rp = c->riph;
      //rm = rp-c->dr;
      //if( i!=0 ) rm = theCells[i-1].riph;
      //double dV = get_dV( rp , rm ); 
      c->dm = c->cons[DDD];//c->prim[RHO]*dV;
   }

   aggregate_mass( theDomain );

}

void grav_src( struct cell * c , double dVdt ){

   double rho = c->prim[RHO];
   double v   = c->prim[VRR];

   double f = get_g( c );

   c->cons[SRR] += rho*f*dVdt;
   if( grav_E_mode == 0 ) c->cons[TAU] += rho*v*f*dVdt;

}

void gravity_addsrc( struct domain * theDomain , double dt ){

   calculate_mass( theDomain );

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int i;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      double rp,rm;
      rp = c->riph;
      rm = rp - c->dr;
      double dV = get_dV( rp , rm );
      grav_src( c , dV*dt );
   }

}


