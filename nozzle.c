
#include "paul.h"

void source_nozz( double * prim , double * cons , double rp , double rm , double t , double dVdt ){

	//double sigma_heat = 7.5e6; //cm
	//double A_heat = 1e15;//1e15; //erg/g/s
	//double r_heat = 4.5e8; //cm

	//double r = .5*(rp+rm);
	//cons[TAU] += A_heat*prim[RHO]*prim[XXX]*exp(-(r-r_heat)*(r-r_heat)/(2.*sigma_heat*sigma_heat))*dVdt;

//////////////////////////////////////////////////
//   double r  = .5*(rp+rm);
//   double R = 0.05;
//   double Vol = 4./3.*M_PI*pow(R,3.);
//   double T = 0.1;
//   double E = 1.0;

//   double Q = 0.0;
//   if( r<R && t<T ) Q = E/T/Vol;

//   cons[TAU] += Q*dVdt;
}


