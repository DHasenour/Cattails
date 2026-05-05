
#include "../paul.h"
#include <mpi.h>

double get_moment_arm( double , double );
double get_dV( double , double );

double get_eint( double * , double * );
double get_entr( double * , double * );
double get_cs( double * , double * );

void calculate_mass( struct domain * );
double get_GMr( struct cell * );

void output( struct domain * theDomain , char * filestart ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Ng = theDomain->Ng;
   int rank = theDomain->rank;
   int size = theDomain->size;

   calculate_mass( theDomain );

   char filename[256];
   sprintf(filename,"%s.dat",filestart);

   if( rank==0 ){
      FILE * pFile = fopen( filename , "w" );
      fprintf(pFile,"#r:1   dr:2   Density:3   Pressure:4   Velocity:5   Alpha:6   X_He:7 - X_Ni:NUM_I+6   M_enc:NUM_I+7   Entropy:NUM_I+8   Temperature:NUM_I+9   cs:NUM_I+10   etot:NUM_I+11\n");
      fclose(pFile);
   }
   MPI_Barrier( MPI_COMM_WORLD );

   int i_min = 0;
   int i_max = Nr;

   if( rank != 0      ) i_min = Ng;
   if( rank != size-1 ) i_max = Nr-Ng;

   double etot , stot , cs;
   int rk;
   for( rk=0 ; rk<size ; ++rk ){
      if( rank==rk ){
         FILE * pFile = fopen( filename , "a" );
         int i,q;
         for( i=i_min ; i<i_max ; ++i ){
            struct cell * c = theCells+i;
            double rp = c->riph;
            double dr = c->dr;
            double rm = rp-dr;
            double r  = .5*(rp+rm);//get_moment_arm( rp , rm );
            fprintf(pFile,"%.14e %.14e ",r,dr);
            for( q=0 ; q<NUM_Q ; ++q ){
               fprintf(pFile,"%.14e ",c->prim[q]);
            }
            fprintf(pFile,"%.14e ",c->miph);

            etot = get_eint( c->prim , &c->T );
            cs = get_cs( c->prim , &c->T );
            stot = get_entr( c->prim , &c->T );
            fprintf(pFile,"%.14e %.14e %.14e %.14e ",stot,c->T,cs,etot);

            /*
            double GMr = get_GMr( c );
            fprintf(pFile,"%.14e ",GMr);

            double ebind = 0.0;
            int q;
            for(q=XXX ; q<NUM_Q ; q++) { ebind += c->prim[RHO]*EBIND[q-XXX]*c->prim[q]; }
            fprintf(pFile,"%.14e ",ebind);
            */

            fprintf(pFile,"\n");
         }
         fclose( pFile );
      }
      MPI_Barrier( MPI_COMM_WORLD );
   }
}
