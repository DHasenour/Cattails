
#include "paul.h"

void reacstep( struct domain * , double );
void onestep( struct domain * , double , double , int , int );

void timestep( struct domain * theDomain , double dt ){
   
   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int i;

   reacstep( theDomain , 0.5*dt );

   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      memcpy( c->RKcons , c->cons , NUM_Q*sizeof(double) );
   }
   onestep( theDomain , 0.0 ,     dt , 1 , 0 );
   onestep( theDomain , 0.5 , 0.5*dt , 0 , 1 );

   reacstep( theDomain , 0.5*dt );

   /////////// 1st Order Strang Splitting //////////////
   //reacstep( theDomain, dt );
   //for( i=0 ; i<Nr ; ++i ){
   //   struct cell * c = theCells+i;
   //   memcpy( c->RKcons , c->cons , NUM_Q*sizeof(double) );
   //}
   //onestep( theDomain , 0.0 ,     dt , 1 , 1 );

   theDomain->t += dt;   
   theDomain->count_steps += 1;

}
