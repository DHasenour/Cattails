
#include "../paul.h"

static double RHO_FLOOR = 0.0;
static double PRE_FLOOR = 0.0;
static int USE_RT = 1;
static double rt_A = 0.0;
static double rt_B = 0.0;
static double rt_C = 0.0;
static double rt_D = 0.0;

void init_helmeos_( char * );
void getcons_( double * , double * , double * , double * , double * , double * , double * , double * );
void getprim_( double * , double * , double * , double * , double * , double * , double * );
void getambient_( double * , double * , double * , double * , double * , double * , double * , double * , double * , double * );
void getpppeos_( double * , double * , double * , double * , double * , double * , double * , double * , double * , double * , double * , double * );
void gettempfromentropy_( double * , double * , double * , double * , double * );

void setHydroParams( struct domain * theDomain ){
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
   USE_RT = theDomain->theParList.rt_flag;
   rt_A = theDomain->theParList.rt_A;
   rt_B = theDomain->theParList.rt_B;
   rt_C = theDomain->theParList.rt_C;
   rt_D = theDomain->theParList.rt_D;
}

void init_eos( char * path2table ){
   //Read in the helm_table (only need to do per execution)
   //Input the directory information
   // "" , "Hydro/Helmeos/" , or "../../Hydro/Helmeos/"
   init_helmeos_( path2table );
}

void azbar( double * X , double * abar , double * zbar ){
   //Routine from Helmeos
   //Given a composition of mass fractions
   //Returns average atomoic number and average charge number
   double a = 0.0;
   double z = 0.0;
   for( int q=0 ; q<NUM_I ; q++ ){
      a += X[q]/aion[q];
      z += X[q]*zion[q]/aion[q];
   }
   *abar = 1.0/a;
   *zbar = z/a;
}

double get_vr( double * prim ){
   //Get radial velocity
   return( prim[VRR] );
}

double get_cs( double * prim , double * T ){
   //Get sound speed given density, pressure, and composition
   double rho = prim[RHO];
   double pre = prim[PPP];
   double x[NUM_I],abar,zbar,etot,stot,temp,cs;
   for( int q=0 ; q<NUM_I ; q++ ){ x[q]=prim[XXX+q]; }
   temp = *T; //a guess
   azbar(x,&abar,&zbar);
   getcons_(&rho,&pre,&abar,&zbar,&etot,&stot,&temp,&cs);
   return cs;
}

double get_pre( double * prim , double * T ){
   //Get pressure given density, temperature, and composition
   double rho = prim[RHO];
   double temp,x[NUM_I],pre,etot,gam,stot,cs,cv,abar,zbar;
   for( int q=0 ; q<NUM_I ; q++ ){ x[q]=prim[XXX+q]; }
   temp = *T;
   azbar(x,&abar,&zbar);
   getambient_(&rho,&temp,&abar,&zbar,&pre,&stot,&gam,&etot,&cs,&cv);
   return pre;
}

double get_pre_from_etot( double * cons , double * T ){
   //Get pressure given density, specific internal energy, and composition
   //cons[density,eint,0,0,mass_frac] : special cons array
   double rho = cons[DDD];
   double etot = cons[TAU];
   double x[NUM_I],abar,zbar,pre,temp,cs;
   for( int q=0 ; q<NUM_I ; q++ ){ x[q]=cons[XXX+q]; }
   temp = *T; //a guess
   azbar(x,&abar,&zbar);
   getprim_(&rho,&etot,&abar,&zbar,&pre,&temp,&cs);
   return pre;
}

double get_temp( double * prim , double * T ){
   //Get temperature given density, pressure, and composition
   double rho = prim[RHO];
   double pre = prim[PPP];
   double x[NUM_I],abar,zbar,etot,stot,temp,cs;
   for( int q=0 ; q<NUM_I ; q++) { x[q]=prim[XXX+q]; }
   temp = *T; //a guess
   azbar(x,&abar,&zbar);
   getcons_(&rho,&pre,&abar,&zbar,&etot,&stot,&temp,&cs);
   return temp;
}

double get_eint( double * prim , double * T ){
   //Get specific internal energy given density, temperature, and composition
   double rho = prim[RHO];
   double temp,x[NUM_I],pre,stot,gam,etot,cs,cv,abar,zbar;
   for( int q=0 ; q<NUM_I ; q++ ){ x[q]=prim[XXX+q]; }
   temp = *T;
   azbar(x,&abar,&zbar);
   getambient_(&rho,&temp,&abar,&zbar,&pre,&stot,&gam,&etot,&cs,&cv);
   return etot;
}

double get_entr( double * prim , double * T ){
   //Get specific entropy given density, temperature, and composition
   double rho = prim[RHO];
   double temp,x[NUM_I],pre,etot,gam,stot,cs,cv,abar,zbar;
   for( int q=0 ; q<NUM_I ; q++ ){ x[q]=prim[XXX+q]; }
   temp = *T;
   azbar(x,&abar,&zbar);
   getambient_(&rho,&temp,&abar,&zbar,&pre,&stot,&gam,&etot,&cs,&cv);
   return stot;
}

double get_temp_cons( double * cons , double * T ){
   //Get temperature given density, specific internal energy, and composition
   //cons[density,eint,0,0,mass_frac] : special cons array
   double rho = cons[DDD];
   double etot = cons[TAU];
   double x[NUM_I],abar,zbar,pre,temp,cs;
   for( int q=0 ; q<NUM_I ; q++) { x[q]=cons[XXX+q]; }
   temp = *T; //a guess
   azbar(x,&abar,&zbar);
   getprim_(&rho,&etot,&abar,&zbar,&pre,&temp,&cs);
   return temp;
}

void get_derivs( double * prim , double * T , double * derivs ){
   //Get dpd, dpt, dsd, and dst given density, temperature, and composition
   //derivs[0] = dpd = change in pressure w.r.t. density
   //derivs[1] = dpt = change in pressure w.r.t. temperature
   //derivs[2] = dsd = change in specific entropy w.r.t. density
   //derivs[3] = dst = change in specific entropy w.r.t. temperature
   double rho = prim[RHO];
   double temp,x[NUM_I],pre,ent,cs,etot,dpd,dpt,dsd,dst,abar,zbar;
   for( int q=0 ; q<NUM_I ; q++ ){ x[q]=prim[XXX+q]; }
   temp = *T;
   azbar(x,&abar,&zbar);
   getpppeos_(&rho,&temp,&abar,&zbar,&pre,&ent,&cs,&etot,&dpd,&dpt,&dsd,&dst);
   
   derivs[0] = dpd;
   derivs[1] = dpt;
   derivs[2] = dsd;
   derivs[3] = dst;
}

void prim2cons( double * prim , double * cons , double GMr , double dV , double * T ){
   int q;
   double rho = prim[RHO];
   double pre = prim[PPP];
   double vr  = prim[VRR];
   double v2 = vr*vr;
   double x[NUM_I],abar,zbar,etot,stot,temp,cs;
   for( q=0 ; q<NUM_I ; q++) { x[q]=prim[XXX+q]; }
   temp = *T; //a guess
   azbar(x,&abar,&zbar);
   getcons_(&rho,&pre,&abar,&zbar,&etot,&stot,&temp,&cs);

   double egrav = -rho*GMr;

   cons[DDD] = rho*dV;
   cons[SRR] = rho*vr*dV;
   cons[TAU] = (.5*rho*v2 + rho*etot + egrav)*dV;
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

   if( rho<RHO_FLOOR ) rho=RHO_FLOOR;
   prim[RHO] = rho;
   prim[VRR] = vr;
   prim[AAA] = cons[AAA]/cons[DDD];

   double x[NUM_I],abar,zbar,pre,etot,temp,cs;
   etot = (E - .5*rho*v2 - egrav)/rho;
   for( q=XXX ; q<NUM_Q ; q++ ){
      prim[q] = cons[q]/cons[DDD];
      etot -= EBIND[q-XXX]*prim[q];
      x[q-XXX]=prim[q];
   }
   temp = *T; //a guess
   azbar(x,&abar,&zbar);
   getprim_(&rho,&etot,&abar,&zbar,&pre,&temp,&cs);
 
   if( pre < PRE_FLOOR*rho ) pre = PRE_FLOOR*rho;
   prim[PPP] = pre;

   //update temperature//
   *T = temp;
}

void getUstar( double * prim , double * Ustar , double Sk , double Ss , double * T ){
   int q;
   double rho = prim[RHO];
   double pre = prim[PPP];
   double vr  = prim[VRR];
   double v2 = vr*vr;
   double x[NUM_I],abar,zbar,etot,stot,temp,cs;
   for( q=0 ; q<NUM_I ; q++ ){ x[q]=prim[XXX+q]; }
   temp = *T; //a guess
   azbar(x,&abar,&zbar);
   getcons_(&rho,&pre,&abar,&zbar,&etot,&stot,&temp,&cs);

   double rhostar = rho*(Sk - vr)/(Sk - Ss);
   double Pstar = pre*(Ss - vr)/(Sk - Ss);
   double Us = .5*rho*v2 + rho*etot;

   Ustar[DDD] = rhostar;
   Ustar[AAA] = prim[AAA]*Ustar[DDD];
   for( q=XXX ; q<NUM_Q ; ++q ){
      Us += rho*EBIND[q-XXX]*prim[q];
      Ustar[q] = prim[q]*Ustar[DDD];
   }
   Us = Us*(Sk - vr)/(Sk - Ss);
   Ustar[SRR] = rhostar*( Ss );
   Ustar[TAU] = Us + rhostar*Ss*(Ss - vr) + Pstar;
}

void flux( double * prim , double * flux , double * T ){
   int q;
   double rho = prim[RHO];
   double pre = prim[PPP];
   double vr  = prim[VRR];
   double v2 = vr*vr;
   double x[NUM_I],abar,zbar,etot,stot,temp,cs;
   for( q=0 ; q<NUM_I ; q++ ){ x[q]=prim[XXX+q]; }
   temp = *T; //a guess
   azbar(x,&abar,&zbar);
   getcons_(&rho,&pre,&abar,&zbar,&etot,&stot,&temp,&cs);
   double rhoe = rho*etot;
 
   flux[DDD] = rho*vr;
   flux[AAA] = flux[DDD]*prim[AAA];
   for( q=XXX ; q<NUM_Q ; ++q ){
      rhoe += rho*EBIND[q-XXX]*prim[q];
      flux[q] = flux[DDD]*prim[q];
   }
   flux[SRR] = rho*vr*vr + pre;
   flux[TAU] = (.5*rho*v2 + rhoe + pre )*vr;
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

