
#include "../paul.h"

static double GAMMA_LAW = 4./3.;

static double RHO_FLOOR = 0.0;
static double PRE_FLOOR = 0.0;
static double grav_G = 0.0;
static int USE_RT = 1;
static double rt_A = 0.0;
static double rt_B = 0.0;
static double rt_C = 0.0;
static double rt_D = 0.0;

void setHydroParams( struct domain * theDomain ){
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
   USE_RT = theDomain->theParList.rt_flag;
   grav_G = theDomain->theParList.grav_G;
   rt_A = theDomain->theParList.rt_A;
   rt_B = theDomain->theParList.rt_B;
   rt_C = theDomain->theParList.rt_C;
   rt_D = theDomain->theParList.rt_D;
}

void init_eos(){
}

double get_vr( double * prim ){
   return( prim[VRR] );
}

double get_cs( double * prim , double * T ){
   //Get sound speed knowing density, pressure, and composition

   return c_light/sqrt(3.);
}

double get_pre( double * prim , double * T ){
   //Get pressure knowing density, temperature, and composition
   double temp = *T;

   return asol*temp*temp*temp*temp/3.;
}

double get_pre_from_etot( double * cons , double * T ){
   //Get pressure knowing density, thermal energy, and composition
   //cons[density,eint,0,0,mass_frac]
   double rho = cons[DDD];
   double rhoe = rho*cons[TAU];

   return (GAMMA_LAW-1.0)*rhoe;
}

double get_temp( double * prim , double * T ){
   //Get temperature knowing density, pressure, and composition
   double pre = prim[PPP];

   return pow(3.*pre/asol,0.25);
}

double get_eint( double * prim , double * T ){
   //Get thermal energy knowing density, temperature, and composition
   double rho = prim[RHO];
   double temp = *T;
   
   return asol*temp*temp*temp*temp/rho;
}

double get_entr( double * prim , double * T ){
   //Get entropy knowing density, temperature, and composition
   double temp = *T;

   return 4./3.*asol*temp*temp*temp;
}

double get_temp_cons( double * cons , double * T ){
   //Get temperature knowing density, thermal energy, and composition
   //cons[density,eint,0,0,mass_frac]
   double rho = cons[DDD];
   double rhoe = rho*cons[TAU];

   return pow(rhoe/asol,0.25);
}

void prim2cons( double * prim , double * cons , double GMr , double dV , double * T ){
   int q;
   double rho = prim[RHO];
   double pre = prim[PPP];
   double vr  = prim[VRR];
   double v2 = vr*vr;

   double rhoe = pre/(GAMMA_LAW-1.);
   double egrav = -rho*GMr;

   cons[DDD] = rho*dV;
   cons[SRR] = rho*vr*dV;
   cons[TAU] = (.5*rho*v2 + rhoe + egrav)*dV;
   cons[AAA] = cons[DDD]*prim[AAA];

   for( q=XXX ; q<NUM_Q ; ++q ){
      cons[TAU] += rho*dV*EBIND[q-XXX]*prim[q];
      cons[q] = cons[DDD]*prim[q];
   }
}

void cons2prim( double * cons , double * prim , double GMr , double dV , double * T ){
   int q;
   double rho = cons[DDD]/dV;
   double Sr  = cons[SRR]/dV;
   double E   = cons[TAU]/dV;
   double egrav = -rho*GMr;
   double vr = Sr/rho;
   double v2 = vr*vr;
   double rhoe = E - .5*rho*v2 - egrav;
   
   if( rho<RHO_FLOOR ) rho=RHO_FLOOR;
   prim[RHO] = rho;
   prim[VRR] = vr;
   prim[AAA] = cons[AAA]/cons[DDD];

   for( q=XXX ; q<NUM_Q ; q++ ){
      prim[q] = cons[q]/cons[DDD];
      rhoe -= rho*EBIND[q-XXX]*prim[q];
   }

   double pre = (GAMMA_LAW-1.)*rhoe;

   if( pre < PRE_FLOOR*rho ) pre = PRE_FLOOR*rho;
   prim[PPP] = pre;

   //update temperature//
   double temp = get_temp( prim , T );
   *T = temp;
}

void getUstar( double * prim , double * Ustar , double Sk , double Ss , double * T ){

   double rho = prim[RHO];
   double vr  = prim[VRR];
   double pre = prim[PPP];
   double v2  = vr*vr;

   double rhoe = pre/(GAMMA_LAW-1.);

   double rhostar = rho*(Sk - vr)/(Sk - Ss);
   double Pstar = pre*(Ss - vr)/(Sk - Ss);
   double Us = rhoe*(Sk - vr)/(Sk - Ss);

   Ustar[DDD] = rhostar;
   Ustar[SRR] = rhostar*( Ss );
   Ustar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss*(Ss - vr) + Pstar;

   int q;
   for( q=XXX ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DDD];
   }

}

void flux( double * prim , double * flux , double * T ){

   double rho = prim[RHO];
   double pre = prim[PPP];
   double vr  = prim[VRR];
   double v2  = vr*vr;

   double rhoe = pre/(GAMMA_LAW-1.);
 
   flux[DDD] = rho*vr;
   flux[SRR] = rho*vr*vr + pre;
   flux[TAU] = (.5*rho*v2 + rhoe + pre )*vr;

   int q;
   for( q=XXX ; q<NUM_Q ; ++q ){
      flux[q] = flux[DDD]*prim[q];
   }

}

void renormalize_comp( double * prim ){
   int q;
   double X_sum = 0.0;
   for( q=XXX ; q<NUM_Q ; ++q ){ X_sum += prim[q]; }
   for( q=XXX ; q<NUM_Q ; ++q ){ prim[q] = prim[q]/X_sum; }
}

void source( double * prim , double * cons , double rp , double rm , double dVdt ){
   double pre = prim[PPP];
   double r  = .5*(rp+rm);
   double r2 = (rp*rp+rm*rm+rp*rm)/3.;
   cons[SRR] += 2.*pre*(r/r2)*dVdt;
}

void source_alpha( double * prim , double * cons , double * grad_prim , double r , double dVdt , double * T ){

   double A = rt_A;//2e-5/1.7; //2e-5;//1e-4;
   double B = rt_B;//1.2;//0.9;
   double D = rt_D;//0.0;

   double alpha = prim[AAA];

   //double pre = prim[PPP];
   double rho = prim[RHO];
   double P1 = grad_prim[PPP];
   double rho1 = grad_prim[RHO];
 
   double g2 = -P1*rho1;
   if( g2 < 0.0 ) g2 = 0.0;
   double cs = get_cs( prim , T );

   cons[AAA] += ( (A+B*alpha)*sqrt(g2) - D*rho*alpha*cs/r )*dVdt;
   if( cons[AAA] < 0. ) cons[AAA] = 0.;

}

double get_eta( double * prim , double * grad_prim , double r , double * T ){

   double C = rt_C;//0.06*1.7;//0.03;

   double cs = get_cs( prim , T );

   double alpha = prim[AAA];
   if( alpha < 0.0 ) alpha = 0.0;

   double u_eddy = cs*sqrt(alpha);
   double lambda = r*sqrt(alpha);
   
   double eta = C*u_eddy*lambda;

   return( eta );

}

void vel( double * prim1 , double * prim2 , double * Sl , double * Sr , double * Ss , double * T1 , double * T2 ){

   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];
   double vn1  = prim1[VRR];

   double cs1 = get_cs( prim1 , T1 );

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];
   double vn2  = prim2[VRR];

   double cs2 = get_cs( prim2 , T2 );

   *Ss = ( P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )/( rho1*(-cs1) - rho2*cs2 );

   *Sr =  cs1 + vn1;
   *Sl = -cs1 + vn1;

   if( *Sr <  cs2 + vn2 ) *Sr =  cs2 + vn2;
   if( *Sl > -cs2 + vn2 ) *Sl = -cs2 + vn2;
   
}

double mindt( double * prim , double w , double r , double g , double dr , double * T ){

   //double rho = prim[RHO];
   //double pre = prim[PPP];
   double vr  = prim[VRR];

//   pre += g*g/8./M_PI/grav_G;

   double cs = get_cs( prim , T );
   double eta = get_eta( prim , NULL , r , T );

   double maxvr = cs + fabs( vr - w );
   double dt = dr/maxvr;
   double dt_eta = dr*dr/eta;
   if( dt > dt_eta && USE_RT ) dt = dt_eta;

   return( dt );

}

