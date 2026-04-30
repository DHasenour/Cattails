
#include "../paul.h"

static double RHO_FLOOR = 0.0;
static double PRE_FLOOR = 0.0;
static double grav_G = 0.0;
static int USE_RT = 1;
static double rt_A = 0.0;
static double rt_B = 0.0;
static double rt_C = 0.0;
static double rt_D = 0.0;

double superCrude_get_Pressure( double , double , double , double );
double superCrude_get_Edens( double , double , double , double );
double superCrude_get_cs2( double , double , double , double );
double superCrude_get_Temp( double , double , double , double );
double superCrude_get_Entropy( double , double , double , double );
double superCrude_get_Temp_fromPre( double , double , double , double );

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

void azbar( double * X , double * abar , double * zbar ){
   int q;
   double a = 0.0;
   double z = 0.0;
   for( q=0 ; q<NUM_I ; q++ ){
      a += X[q]/aion[q];
      z += X[q]*zion[q]/aion[q];
   }
   *abar = 1.0/a;
   *zbar = z/a;
}

double get_vr( double * prim ){
   return( prim[VRR] );
}

double get_cs( double * prim , double * T ){
   //Get sound speed knowing density, pressure, and composition
   int q;
   double rho = prim[RHO];
   double pre = prim[PPP];
   double x[NUM_I],abar,zbar,temp,cs2;
   for( q=0 ; q<NUM_I ; q++) { x[q]=prim[XXX+q]; }
   azbar(x,&abar,&zbar);
   temp = superCrude_get_Temp_fromPre(rho,pre,abar,zbar);
   cs2 = superCrude_get_cs2(rho,temp,abar,zbar);
   return sqrt(cs2);
}

double get_pre( double * prim , double * T ){
   //Get pressure knowing density, temperature, and composition
   int q;
   double rho = prim[RHO];
   double temp,x[NUM_I],pre,abar,zbar;
   for( q=0 ; q<NUM_I ; q++ ){ x[q]=prim[XXX+q]; }
   temp = *T;
   azbar(x,&abar,&zbar);
   pre = superCrude_get_Pressure(rho,temp,abar,zbar);
   return pre;
}

double get_pre_from_etot( double * cons , double * T ){
   //Get pressure knowing density, thermal energy, and composition
   //cons[density,eint,0,0,mass_frac]
   int q;
   double rho = cons[DDD];
   double etot = cons[TAU];
   double x[NUM_I],abar,zbar,pre,temp;
   for( q=0 ; q<NUM_I ; q++) { x[q]=cons[XXX+q]; }
   azbar(x,&abar,&zbar);
   temp = superCrude_get_Temp(rho,rho*etot,abar,zbar);
   pre = superCrude_get_Pressure(rho,temp,abar,zbar);
   return pre;
}

double get_temp( double * prim , double * T ){
   //Get temperature knowing density, pressure, and composition
   int q;
   double rho = prim[RHO];
   double pre = prim[PPP];
   double x[NUM_I],abar,zbar,temp;
   for( q=0 ; q<NUM_I ; q++) { x[q]=prim[XXX+q]; }
   azbar(x,&abar,&zbar);
   temp = superCrude_get_Temp_fromPre(rho,pre,abar,zbar);
   return temp;
}

double get_eint( double * prim , double * T ){
   //Get thermal energy knowing density, temperature, and composition
   int q;
   double rho = prim[RHO];
   double temp,x[NUM_I],rhoe,abar,zbar;
   for( q=0 ; q<NUM_I ; q++ ){ x[q]=prim[XXX+q]; }
   temp = *T;
   azbar(x,&abar,&zbar);
   rhoe = superCrude_get_Edens(rho,temp,abar,zbar);
   return rhoe / rho;
}

double get_entr( double * prim , double * T ){
   //Get entropy knowing density, temperature, and composition
   int q;
   double rho = prim[RHO];
   double temp,x[NUM_I],stot,abar,zbar;
   for( q=0 ; q<NUM_I ; q++ ){ x[q]=prim[XXX+q]; }
   temp = *T;
   azbar(x,&abar,&zbar);
   stot = superCrude_get_Entropy(rho,temp,abar,zbar);
   return stot;
}

double get_temp_cons( double * cons , double * T ){
   //Get temperature knowing density, thermal energy, and composition
   //cons[density,eint,0,0,mass_frac]
   int q;
   double rho = cons[DDD];
   double etot = cons[TAU];
   double x[NUM_I],abar,zbar,temp;
   for( q=0 ; q<NUM_I ; q++) { x[q]=cons[XXX+q]; }
   azbar(x,&abar,&zbar);
   temp = superCrude_get_Temp(rho,rho*etot,abar,zbar);
   return temp;
}

void prim2cons( double * prim , double * cons , double GMr , double dV , double * T ){
   int q;
   double rho = prim[RHO];
   double pre = prim[PPP];
   double vr  = prim[VRR];
   double v2 = vr*vr;
   double x[NUM_I],abar,zbar,rhoe,temp;
   for( q=0 ; q<NUM_I ; q++) { x[q]=prim[XXX+q]; }
   azbar(x,&abar,&zbar);
   temp = superCrude_get_Temp_fromPre(rho,pre,abar,zbar);
   rhoe = superCrude_get_Edens(rho,temp,abar,zbar);

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

   if( rho<RHO_FLOOR ) rho=RHO_FLOOR;
   prim[RHO] = rho;
   prim[VRR] = vr;
   prim[AAA] = cons[AAA]/cons[DDD];

   double x[NUM_I],abar,zbar,pre,etot,temp;
   etot = (E - .5*rho*v2 - egrav)/rho;
   for( q=XXX ; q<NUM_Q ; q++ ){
      prim[q] = cons[q]/cons[DDD];
      etot -= EBIND[q-XXX]*prim[q];
      x[q-XXX]=prim[q];
   }
   azbar(x,&abar,&zbar);
   temp = superCrude_get_Temp(rho,rho*etot,abar,zbar);
   pre = superCrude_get_Pressure(rho,temp,abar,zbar);
 
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
   double x[NUM_I],abar,zbar,rhoe,temp;
   for( q=0 ; q<NUM_I ; q++ ){ x[q]=prim[XXX+q]; }
   azbar(x,&abar,&zbar);
   temp = superCrude_get_Temp_fromPre(rho,pre,abar,zbar);
   rhoe = superCrude_get_Edens(rho,temp,abar,zbar);

   double rhostar = rho*(Sk - vr)/(Sk - Ss);
   double Pstar = pre*(Ss - vr)/(Sk - Ss);
   double Us = .5*rho*v2 + rhoe;

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
   double x[NUM_I],abar,zbar,rhoe,temp;
   for( q=0 ; q<NUM_I ; q++ ){ x[q]=prim[XXX+q]; }
   azbar(x,&abar,&zbar);
   temp = superCrude_get_Temp_fromPre(rho,pre,abar,zbar);
   rhoe = superCrude_get_Edens(rho,temp,abar,zbar);
 
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

//////////////////////////////////////////////////////////////////////////
//                           SUPERCRUDE EOS
//////////////////////////////////////////////////////////////////////////

double superCrude_TOL = 1e-8;

//double kB = 1.380649e-16; //in paul.h (k_B)
//double me = 9.10938e-28; //in paul.h (m_e)
//double mp = 1.67262192e-24; //in paul.h (m_p)
//double c  = 2.998e10; //in paul.h (c_light)
//double hbar = 1.05457266e-27; //in paul.h (hbar)
//double aR = 7.5646e-15; //in paul.h (asol)

int include_deg_factor = 1;

double Integral(double x){
   double I = ( x*sqrt(1.+x*x)*(1.+2.*x*x) - asinh(x) )/8. - x*x*x/3.;
   double x0 = 3e-3;
   double f = exp(-x0/x);
   I = I*f + (1.-f)*0.1*pow(x,5.);  //Replace with Taylor expansion in limit x->0
   return( I );
}

double Integral2( double x ){
   double I = ( x*sqrt(1.+x*x)*(1.+2.*x*x) - asinh(x) )/8. - x*x*x/3.*sqrt(1.+x*x);
   double x0 = 3e-3;
   double f = exp(-x0/x);
   I = I*f - (1.-f)*x*x*x*x*x*2./30.;
   return( I );
}

double get_fermi( double rho , double T , double A , double Z , int mode ){

   double ne;
   if( mode == 0 ) ne = Z*rho*n_A/A; else ne = rho*n_A/A;
   double pF = pow(3.*M_PI*M_PI*ne,1./3.)*hbar;
   double x = pF/m_e/c_light;
   double I = Integral(x);
 
   double eF = m_e*c_light*c_light*( sqrt(1.+x*x)-1. );
   
   return(eF);

}

double P_deg( double rho , double A , double Z ){

   double ne = Z*rho*n_A/A;
   double pF = pow(3.*M_PI*M_PI*ne,1./3.)*hbar;
   double x = pF/m_e/c_light;
   double I = -Integral2(x);

   double e0 = c_light*pow(m_e*c_light,4.)/M_PI/M_PI/pow(hbar,3.);

   return( e0*I );
}

double e_deg( double rho , double A , double Z ){

   double ne = Z*rho*n_A/A;
   double pF = pow(3.*M_PI*M_PI*ne,1./3.)*hbar;
   double x = pF/m_e/c_light;
   double I = Integral(x);

   double e0 = c_light*pow(m_e*c_light,4.)/M_PI/M_PI/pow(hbar,3.);

   return(e0*I);
}

double superCrude_get_Pressure( double rho , double T , double A , double Z ){

   double eF = get_fermi( rho , T , A , Z , 0 );
   double T0 = 6./M_PI/M_PI*eF/k_B;

   double elec_plus_ions = 1.+Z;
if( include_deg_factor )   elec_plus_ions = 1. + Z*(T/(T+T0));
   double n = elec_plus_ions*rho*n_A/A;

   double Pgas = n*k_B*T;
   double Prad = asol*T*T*T*T/3.;
   double Pdeg = P_deg( rho , A , Z );

//   double b = m_e*c_light*c_light/k_B/T;
//   double Ib = sqrt(M_PI/2.)/pow(b,1.5) + 2./b/b/b;
//   double npos = 2.*pow(m_e*c_light/hbar,3.)/M_PI/M_PI*Ib*exp(-b);
//   double Ppos = npos*k_B*T;
   
   return( Pgas + Prad + Pdeg );//+ Ppos );

}

double superCrude_get_Edens( double rho , double T , double A , double Z ){

   double eF = get_fermi( rho , T , A , Z , 0 );
   double T0 = 6./M_PI/M_PI*eF/k_B;

   double elec_plus_ions = 1.+Z;
if( include_deg_factor )   elec_plus_ions = 1. + Z*T/(T+T0); //*( 1. - exp(-T/T0) );
   double n = elec_plus_ions*rho*n_A/A;

   double egas = 1.5*n*k_B*T;
   double erad = asol*T*T*T*T;
   double edeg = e_deg( rho , A , Z );

//   double b = m_e*c*c/k_B/T;
//   double Ib = sqrt(M_PI/2.)/pow(b,1.5) + 2./b/b/b;
//   double npos = 2.*pow(m_e*c/hbar,3.)/M_PI/M_PI*Ib*exp(-b);
//   double epos = (3. + b)*npos*k_B*T;

   return( egas + erad + edeg );//+ epos );

}

double superCrude_get_cs2( double rho , double T , double A , double Z ){
   
   double eF = get_fermi( rho , T , A , Z , 0 );
   double T0 = 6./M_PI/M_PI*eF/k_B;
   
   double elec_plus_ions = 1.+Z;
if( include_deg_factor )   elec_plus_ions = 1. + Z*(T/(T+T0));
   double n = elec_plus_ions*rho*n_A/A;

//   double gamPgas = 5./3.*n*k_B*T;
//   double gamPrad = 4./3.*asol*T*T*T*T/3.;

   double gamPgasrad = n*k_B*T + pow(n*k_B*T + 4./3.*asol*T*T*T*T ,2.)/( 1.5*n*k_B*T + 4.*asol*T*T*T*T );

   double ne = Z*rho*n_A/A;
   double pF = hbar*pow( 3.*M_PI*M_PI*ne , 1./3. );
   double gamPdeg = 1./3.*ne*pF*pF*c_light/sqrt( m_e*m_e*c_light*c_light + pF*pF );

   double P   = superCrude_get_Pressure( rho , T , A , Z );
   double eps = superCrude_get_Edens( rho , T , A , Z );

   double rhoh = rho + (P+eps)/c_light/c_light;

   double cs2 = (gamPgasrad + gamPdeg)/rhoh;

   return( cs2 );

}

double E_gasrad( double rho , double T , double A , double Z ){

   double eF = get_fermi( rho , T , A , Z , 0 );
   double T0 = 6./M_PI/M_PI*eF/k_B;

   double elec_plus_ions = 1.+Z;
   if( include_deg_factor ) elec_plus_ions = 1. + Z*(T/(T+T0));
   double n = elec_plus_ions*rho*n_A/A;

   double egas = 1.5*n*k_B*T;
   double erad = asol*T*T*T*T;

   return( egas + erad );

}

double E_deriv( double rho , double T , double A , double Z ){

   double eF = get_fermi( rho , T , A , Z , 0 );
   double T0 = 6./M_PI/M_PI*eF/k_B;

   double elec_plus_ions = 1.+Z;
   if( include_deg_factor ) elec_plus_ions = 1. + Z*(T/(T+T0));
   double n = elec_plus_ions*rho*n_A/A;
   double ne = Z*rho*n_A/A;

   double e1gas = 1.5*n*k_B;
   if( include_deg_factor ) e1gas += 1.5*ne*k_B*T*T0/pow(T+T0,2.);
   double e1rad = 4.*asol*T*T*T;

   return( e1gas + e1rad );

}

double get_Tguess( double rho , double eps , double A , double Z ){

   double elec_plus_ions = Z+1.;
   double n = elec_plus_ions*rho*n_A/A;

   double T1 = eps/(1.5*n*k_B);
   double T2 = pow(eps/asol,0.25);

   double Tguess = 1./sqrt( 1./T1/T1 + 1./T2/T2 );

   return(Tguess);

}

double superCrude_get_Temp( double rho , double eps , double A , double Z ){

   double ed = e_deg( rho , A , Z );
   eps -= ed;
   if( eps<0.0 ) return(0.0);

   double Tguess = get_Tguess( rho , eps , A , Z );
   //Tguess = 1e10;
   double eps_guess = E_gasrad( rho , Tguess , A , Z );
   double de = eps_guess-eps;
//   int n = 0;

   while( fabs(de/eps) > superCrude_TOL ){
      double e1 = E_deriv( rho , Tguess , A , Z );
      double dT = -de/e1;
      Tguess += dT;
      eps_guess = E_gasrad( rho , Tguess , A , Z );
      de = eps_guess-eps;
//      ++n;
   }
//   printf("%e, %d ",Tguess,n);

   return( Tguess );
}

double superCrude_get_Cv( double rho , double T , double A , double Z ){

   double e1 = E_deriv( rho , T , A , Z );
   return( e1/rho );

}

double superCrude_get_Entropy( double rho , double T , double A , double Z ){
   
   double eF = get_fermi( rho , T , A , Z , 0 );
   double T0 = 6./M_PI/M_PI*eF/k_B;

   double ne = Z*rho*n_A/A;
   double he = 2.5;
   if( include_deg_factor ) he *= T/(T0+T); //(1.-exp(-T/T0));
//   if( include_deg_factor ) ne *= ( 1. - exp(-T/T0) );
   double ni = rho*n_A/A;

   double Srad = 4./3.*asol*T*T*T;

   double Sele = ne*k_B*( he  + 1.5*log( 1. + M_PI*M_PI/18.*k_B*T/eF ) );
   double Sion = ni*k_B*( 1.5*log( 1. + exp(2.5)*(A/n_A)*k_B*T/2./M_PI/hbar/hbar/pow(ni,2./3.) ) );

   double s = ( Srad + Sele + Sion )/rho;
//   if( s<0.0 ) s=0.0;

   return( s );

}

double superCrude_get_Temp_fromPre( double rho , double pre , double A , double Z ){

   return( 0.0 );
}
