
#include "../../paul.h"

//functions from euler-'EOS'.c
void init_eos( int );
double get_pre( double * , double * );
double get_entr( double * , double * );
void get_derivs( double * , double * , double * );

double getPPPbalance( double , double , double , double * , double );
void getPPPandENTbalance( double , double , double , double * , double , double * , double * );
void getdensities( double , double , double * , double * );
void getmasses( double , double , double * , double * , int );

///////////////////////////////////////////////////////////
///////////////////////   MAIN   //////////////////////////
///////////////////////////////////////////////////////////

int main( int argc, char *argv[] ){
   init_eos(2);

   //mode: 0 = WD finder ; 1 = WD init 
   int mode = 0;

   //Desired White Dwarf properties
   double mass_total = 1.00 * Msun;
   double mass_shell = 0.050 * Msun;

   //Initial Guesses
   double den_core = 3.33374669687980e+07;
   double den_He = 1.24201160233603e+06;

   if( mode == 0 ){
      getdensities( mass_total , mass_shell , &den_core , &den_He );
      //if converged... write out WD
      mode = 1;
   }

   if( mode == 1 ){
      getmasses( den_core , den_He , &mass_total , &mass_shell , mode );
   }

   return(0);
}

///////////////////////////////////////////////////////////
/////////////////////// END MAIN //////////////////////////
///////////////////////////////////////////////////////////

void getmasses( double den_core, double den_shell , double * mas_total , double * mas_shell , int flag ){

   int q; 
   ////////// Integration params ///////////
   double temp_core = 1e7; //temperature of the isothermal core
   double Nr = 1048576; //total number of radial zones
   double Rmin = 0.0; //inner boundary of WD
   double Rmax = 2e11; //outer boundary of WD
   double R0 = 1e9; //transition point lin to log
   double f = 0.5; //fraction of zones in linear region
   double n = 2.0; //smoothness param
   double den_He = den_shell;
   double delta = 5e6; //width of tanh transition
   double temp_base = 2e8;//1.75e8; //temperature at the base of the shell
   double NUM_delta = 4; //tanh is 2*NUM_delta wide. ie +/- number of widths
   double temp_csmtran = 1e7; //temperature to stop isentropic shell
   double den_csmtran = 0; //idk
   double den_cut = 1e-3; //idk
   /////////////////////////////////////////

   char out_name[256] = "initial.dat";
   FILE * pFile = fopen(out_name,"w"); 

   /////////////////////////////////////////
   // Other initial params // 
   double den = den_core;
   double temp = temp_core;
   double comp[NUM_I],comp_core[NUM_I],comp_shell[NUM_I];
   int index_He4 = -1, index_C12 = -1, index_O16 = -1;
   for( q=0 ; q<NUM_I ; q++ ){ 
      comp_core[q] = 0.0; 
      comp_shell[q] = 0.0;
      if( aion[q] == 4.0 && zion[q] == 2.0 ){ index_He4 = q; }
      if( aion[q] == 12.0 && zion[q] == 6.0 ){ index_C12 = q; }
      if( aion[q] == 16.0 && zion[q] == 8.0 ){ index_O16 = q; }
   }
   if( index_He4 < 0 || index_C12 < 0 || index_O16 < 0 ){
      printf("Error: He4, C12, or O16 is not in network\n");
   }
   
   comp_core[index_C12] = 0.5;
   comp_core[index_O16] = 0.5;

   comp_shell[index_He4] = 1.0;
   /////////////////////////////////////////

   double dx = 1./Nr;
   double rp, rm, dr, Rlin, Rlog, scal;
   double x0 = 0.0;
   double xm = x0;
   double xp = x0 + dx;
   scal = pow(pow((R0-Rmin)/f,n)+pow(Rmax-R0*pow(Rmax/R0,-f/(1.-f)),n),1./n);
   //Outer radius// 
   Rlin = xp*(R0-Rmin)/f;
   Rlog = R0*pow(Rmax/R0,(xp-f)/(1.-f)) - R0*pow(Rmax/R0,-f/(1.-f));
   rp = (Rmax-Rmin)*pow(pow(Rlin,n)+pow(Rlog,n),1./n)/scal + Rmin;
   //Inner radius// 
   Rlin = xm*(R0-Rmin)/f;
   Rlog = R0*pow(Rmax/R0,(xm-f)/(1.-f)) - R0*pow(Rmax/R0,-f/(1.-f));
   rm = (Rmax-Rmin)*pow(pow(Rlin,n)+pow(Rlog,n),1./n)/scal + Rmin;
   dr = rp-rm;

   double mas = 0.;
   double masCO = 0., masHe = 0.;
   double r_imh = 0., r_iph=dr;
   double g_imh; 

   double vel = 0.0;
   double X = 0.0;  //, X2 = 0.0;
   double rr, r3, r2, r_core = 2*Rmax, r_trans = 2*Rmax, r_shell = 2*Rmax, r_trans2 = 2*Rmax;
   double S_He = 0.0;
   double ptot , gam , etot , cs , cv;
   double abar , zbar;
   double prim[NUM_Q];
   double pre;

   // STEP 1: get core pressure //
   prim[RHO] = den;
   for( q=0 ; q<NUM_I ; q++ ){ prim[XXX+q]=comp_core[q]; }
   pre = get_pre( prim , &temp );
   // STEP 1b: update/output core values 
   double den_im1, pre_im1, delta_r, dV;

   delta_r = r_iph-r_imh;
   dV = (4.*M_PI*((r_iph*r_iph+r_imh*r_imh+r_iph*r_imh)/3.)*delta_r);
   mas += den*dV; 
   masCO += den*dV*(1.0-X);
   masHe += den*dV*X;

   //WRITE OUT//
   if( flag == 1 ){
      //r3 = (r_iph*r_iph*r_iph + r_iph*r_iph*r_imh + r_iph*r_imh*r_imh + r_imh*r_imh*r_imh)/4.;
      //r2 = (r_iph*r_iph + r_iph*r_imh + r_imh*r_imh)/3.;
      //rr = r3/r2;
      rr = 0.5*(r_iph+r_imh);
      fprintf(pFile,"%.14e %.14e %.14e %.14e %.14e\n",rr,den,pre,vel,X);
   }

   // STEP 2: Integrate outward //
   for( int i=1 ; i<Nr ; ++i ){

      xm = x0 + ((double)i)*dx;
      xp = x0 + ((double)i+1.)*dx;
      //Outer radius// 
      Rlin = xp*(R0-Rmin)/f;
      Rlog = R0*pow(Rmax/R0,(xp-f)/(1.-f)) - R0*pow(Rmax/R0,-f/(1.-f));
      rp = (Rmax-Rmin)*pow(pow(Rlin,n)+pow(Rlog,n),1./n)/scal + Rmin;
      //Inner radius// 
      Rlin = xm*(R0-Rmin)/f;
      Rlog = R0*pow(Rmax/R0,(xm-f)/(1.-f)) - R0*pow(Rmax/R0,-f/(1.-f));
      rm = (Rmax-Rmin)*pow(pow(Rlin,n)+pow(Rlog,n),1./n)/scal + Rmin;
      
      r_imh = rm;
      r_iph = rp;
      if( r_iph > Rmax ) {r_iph = Rmax;} // special last step
      delta_r = r_iph-r_imh;

      den_im1 = den;
      pre_im1 = pre;

      if( r_iph < r_core ) {
         ////////////////////////  CORE   //////////////////////////////
         for ( q=0 ; q<NUM_I ; q++ ){ comp[q]=comp_core[q]; }
         temp = temp_core;
         g_imh = -grav_G*mas/(r_imh*r_imh);
         den = getPPPbalance( pre_im1 , den_im1 , temp , comp , delta_r*g_imh );
         
         prim[RHO] = den;
         for( q=0 ; q<NUM_I ; q++ ){ prim[XXX+q]=comp[q]; }
         pre = get_pre( prim , &temp );

         if( den < den_He ){
            r_core = r_iph;
            if( flag == 1 ){ printf("CORE  ENDS --- r_core: %e\n",r_core); }
         }
      }
      else if( r_iph < r_trans ){
         //////////////////////  CORE TO SHELL TRANSITION  /////////////////////////////
         X = 0.5*(1.0+tanh((r_iph-r_core-NUM_delta*delta)/delta));
         if( r_iph-r_core > 2*NUM_delta*delta ){ X = 1.0; }
         else { X = (2*X-1.0)/(2*tanh(NUM_delta)) + 0.5; }

         for ( q=0 ; q<NUM_I ; q++ ){comp[q]=comp_core[q] + X*(comp_shell[q] - comp_core[q]);}
         temp = temp_core + X*(temp_base - temp_core);

         g_imh = -grav_G*mas/(r_imh*r_imh);
         den = getPPPbalance( pre_im1 , den_im1 , temp , comp , delta_r*g_imh );

         prim[RHO] = den;
         for( q=0 ; q<NUM_I ; q++ ){ prim[XXX+q]=comp[q]; }
         pre = get_pre( prim , &temp );

         if( r_iph-r_core > 2*NUM_delta*delta ){
            r_trans = r_iph;
            S_He = get_entr( prim , &temp );
            if( flag == 1 ){ printf("TANH  ENDS --- r_trans , S_He: %e %e\n",r_trans,S_He); }
         }
      }
      else if( r_iph < r_shell ){
         //////////////////////  HELIUM SHELL  /////////////////////////////
         for( q=0 ; q<NUM_I ; q++ ){comp[q]=comp_shell[q];}

         g_imh = -grav_G*mas/(r_imh*r_imh);
         getPPPandENTbalance( pre_im1 , den_im1 , S_He , comp , delta_r*g_imh , &den , &temp );
         
         prim[RHO] = den;
         for( q=0 ; q<NUM_I ; q++ ){ prim[XXX+q]=comp[q]; }
         pre = get_pre( prim , &temp );

         if( temp < temp_csmtran ){
            r_shell = r_iph;
            temp_csmtran = temp;
            den_csmtran = den;
            if( flag == 1 ){ printf("ISENT ENDS --- r_shell, den, temp, mass, pre, T_delta: %e %e %e %e %e\n",r_shell,den,temp,mas,pre); }
         }
      }
      else if( r_iph < r_trans2 ){
         //////////////////////  SHELL TO CSM TRANSITION  /////////////////////////////
         for ( q=0 ; q<NUM_I ; q++ ){comp[q]=comp_shell[q];}
         temp = temp_csmtran*pow(r_shell/r_iph,0.75);
         //temp = temp_csmtran;

         g_imh = -grav_G*mas/(r_imh*r_imh);
         den = getPPPbalance( pre_im1 , den_im1 , temp , comp , delta_r*g_imh );
         
         prim[RHO] = den;
         for( q=0 ; q<NUM_I ; q++ ){ prim[XXX+q]=comp[q]; }
         pre = get_pre( prim , &temp );

         if( den < den_cut ){
            //r_trans2 = r_iph;
            //den_cut = den;
            //printf("EXP   ENDS --- r_trans2, den: %e %e\n",r_trans2,den);
         }
      }
      else{
         //////////////////////  CSM  /////////////////////////////
         for ( q=0 ; q<NUM_I ; q++ ){comp[q]=comp_shell[q];}
         den = den_cut;
         pre = pre_im1;
      }
      

      // STEP 2d: accumulate mass
      dV = (4.*M_PI*((r_iph*r_iph+r_imh*r_imh+r_iph*r_imh)/3.)*delta_r);
      mas += den*dV; 
      masCO += den*dV*(1.0-X);
      masHe += den*dV*X;

      //WRITE OUT//
      if( flag == 1 ){
         //r3 = (r_iph*r_iph*r_iph + r_iph*r_iph*r_imh + r_iph*r_imh*r_imh + r_imh*r_imh*r_imh)/4.;
         //r2 = (r_iph*r_iph + r_iph*r_imh + r_imh*r_imh)/3.;
         //rr = r3/r2;
         rr = 0.5*(r_iph+r_imh);
         fprintf(pFile,"%.14e %.14e %.14e %.14e %.14e\n",rr,den,pre,vel,X);
      }
   }

   fclose(pFile);
     
   *mas_total = mas;
   *mas_shell = masHe;
}

void getdensities( double mass_total , double mass_shell , double * den_core , double * den_He ){
   //////////////////////////////////////
   //      Commence Chaos
   //////////////////////////////////////
   int i; 
   double eostol = 1e-6;
   double fpmin = 1e-14; 
   double ep = 0.0001;
   double eoswrk01, eoswrk02, eoswrk03, denCorenew, denShellnew;
   double denCore = *den_core; //initial guess
   double denShell = *den_He; //initial guess
   double masTotal,masTotal1,masTotal2,masTotal3,masTotal4;
   double masShell,masShell1,masShell2,masShell3,masShell4;

   // Manually do inital step //
   eoswrk01 = 0.0;
   eoswrk02 = 0.0;
   eoswrk03 = 0.0;

   getmasses( denCore , denShell , &masTotal , &masShell , 0 );
   getmasses( (1+ep)*denCore , denShell , &masTotal1 , &masShell1 , 0 );
   getmasses( (1-ep)*denCore , denShell , &masTotal2 , &masShell2 , 0 );
   getmasses( denCore , (1+ep)*denShell , &masTotal3 , &masShell3 , 0 );
   getmasses( denCore , (1-ep)*denShell , &masTotal4 , &masShell4 , 0 );

   double f1 = masTotal - mass_total;
   double f2 = masShell - mass_shell;
   double df11 = (masTotal1 - masTotal2)/(2.*ep*denCore);
   double df12 = (masTotal3 - masTotal4)/(2.*ep*denShell);
   double df21 = (masShell1 - masShell2)/(2.*ep*denCore);
   double df22 = (masShell3 - masShell4)/(2.*ep*denShell);
   double det = df11*df22-df12*df21;
   eoswrk02 = (df22*f1-df12*f2)/det;
   eoswrk03 = (df11*f2-df21*f1)/det;

   denCorenew = fmin(fmax(.5*denCore,denCore - eoswrk02),2.*denCore);
   denShellnew = fmin(fmax(.5*denShell,denShell - eoswrk03),2.*denShell);

   eoswrk01 = fabs((denCorenew - denCore)/denCore) + fabs((denShellnew - denShell)/denShell);
   //Helmeos Bounds//
   denCore = fmin(1e14,fmax(denCorenew,1e-11));
   denShell = fmin(1e14,fmax(denShellnew,1e-11));

   int converged = 0;
   // Loop over the remaining steps //
   for( i=2 ; i<101 ; i++ ){
      if( eoswrk01 < eostol || (fabs(eoswrk02)+fabs(eoswrk03)) < fpmin ) {
         //converged!!!
         converged = 1;
         break;
      }

      getmasses( denCore , denShell , &masTotal , &masShell , 0 );
      printf("%d: %.14e %.14e , %.8e %.8e , %e\n",i,denCore,denShell,masTotal/Msun,masShell/Msun,eoswrk01);
      
      getmasses( (1+ep)*denCore , denShell , &masTotal1 , &masShell1 , 0 );
      getmasses( (1-ep)*denCore , denShell , &masTotal2 , &masShell2 , 0 );
      getmasses( denCore , (1+ep)*denShell , &masTotal3 , &masShell3 , 0 );
      getmasses( denCore , (1-ep)*denShell , &masTotal4 , &masShell4 , 0 );

      f1 = masTotal - mass_total;
      f2 = masShell - mass_shell;
      df11 = (masTotal1 - masTotal2)/(2.*ep*denCore);
      df12 = (masTotal3 - masTotal4)/(2.*ep*denShell);
      df21 = (masShell1 - masShell2)/(2.*ep*denCore);
      df22 = (masShell3 - masShell4)/(2.*ep*denShell);
      det = df11*df22-df12*df21;
      eoswrk02 = (df22*f1-df12*f2)/det;
      eoswrk03 = (df11*f2-df21*f1)/det;

      denCorenew = fmin(fmax(.5*denCore,denCore - eoswrk02),2.*denCore);
      denShellnew = fmin(fmax(.5*denShell,denShell - eoswrk03),2.*denShell);

      eoswrk01 = fabs((denCorenew - denCore)/denCore) + fabs((denShellnew - denShell)/denShell);
      denCore = fmin(1e14,fmax(denCorenew,1e-11));
      denShell = fmin(1e14,fmax(denShellnew,1e-11));
   }
   if( converged == 0 ){ printf("ERROR:failed to converge\n"); }
   else{ 
      printf("Final: %.14e %.14e\n",denCore,denShell); 
      *den_core = denCore;
      *den_He = denShell;
   }
}

double getPPPbalance( double pre_im1 , double den_im1 , double temp , double comp[] , double drgimh ) {

   int i; 
   double eostol = 1e-8;
   double fpmin = 1e-14; 
   double eoswrk01, eoswrk02, preHSE, dennew;
   double den = den_im1; //initial guess
   double preEOS, dpd;
   double prim[NUM_Q], derivs[4]; //derivs = [dpd,dpt,dsd,dst]

   // Manually do inital step //
   eoswrk01 = 0.0;
   eoswrk02 = 0.0;
   preHSE = pre_im1 + .5*(den+den_im1)*drgimh; //needs to be HSE pres

   prim[RHO] = den;
   for( i=0 ; i<NUM_I ; i++ ){ prim[XXX+i]=comp[i]; }
   preEOS = get_pre( prim , &temp );
   get_derivs(prim , &temp , derivs);
   dpd = derivs[0];

   double f = preEOS - preHSE;
   double df = dpd - .5*drgimh;
   eoswrk02 = f/df;

   dennew = fmin(fmax(.5*den,den - eoswrk02),2.*den);
   eoswrk01 = fabs((dennew - den)/den);
   den = fmin(1e14,fmax(dennew,1e-11));

   // Loop over the remaining steps //
   for( i=2 ; i<41 ; i++ ){
      if( eoswrk01 < eostol || fabs(eoswrk02) < fpmin ) {
         //converged!!!
         return den;
      }

      preHSE = pre_im1 + .5*(den+den_im1)*drgimh; // HSE pressure
      prim[RHO] = den;
      preEOS = get_pre( prim , &temp );
      get_derivs(prim , &temp , derivs);
      dpd = derivs[0];

      f = preEOS - preHSE;
      df = dpd - .5*drgimh;
      eoswrk02 = f/df;

      dennew = fmin(fmax(.5*den,den - eoswrk02),2.*den);
      eoswrk01 = fabs((dennew - den)/den);
      den = fmin(1e14,fmax(dennew,1e-11));
   }
   // if here... not converged :(
   printf("ERROR: Pressure failed to balance\n");
   return -1.;
}

void getPPPandENTbalance( double pre_im1 , double den_im1 , double S_He , double comp[] , double drgimh , double *density , double *temp ){

   int i; 
   double eostol = 1e-8;
   double fpmin = 1e-14; 
   double eoswrk01, eoswrk02, eoswrk03, preHSE, dennew, tempnew;
   double den = den_im1; //initial guess
   double T = *temp; //initial guess
   double preEOS, entEOS, dpd, dpt, dsd, dst;
   double prim[NUM_Q], derivs[4]; //derivs = [dpd,dpt,dsd,dst]

   // Manually do inital step //
   eoswrk01 = 0.0;
   eoswrk02 = 0.0;
   preHSE = pre_im1 + .5*(den+den_im1)*drgimh; //needs to be HSE pres

   prim[RHO] = den;
   for( i=0 ; i<NUM_I ; i++ ){ prim[XXX+i]=comp[i]; }
   preEOS = get_pre( prim , &T );
   entEOS = get_entr( prim , &T );
   get_derivs(prim , &T , derivs);
   dpd = derivs[0];
   dpt = derivs[1];
   dsd = derivs[2];
   dst = derivs[3];

   double f1 = preEOS - preHSE;
   double f2 = entEOS - S_He;
   double df11 = dpd - .5*drgimh;
   double df12 = dpt;
   double df21 = dsd;
   double df22 = dst;
   double det = df11*df22-df12*df21;
   eoswrk02 = (df22*f1-df12*f2)/det;
   eoswrk03 = (df11*f2-df21*f1)/det;

   dennew = fmin(fmax(.5*den,den - eoswrk02),2.*den);
   tempnew = fmin(fmax(.5*T,T - eoswrk03),2.*T);

   eoswrk01 = fabs((dennew - den)/den) + fabs((tempnew - T)/T);
   den = fmin(1e14,fmax(dennew,1e-11));
   T = fmin(1e13,fmax(tempnew,1e3));

   // Loop over the remaining steps //
   for( i=2 ; i<101 ; i++ ){
      if( eoswrk01 < eostol || (fabs(eoswrk02)+fabs(eoswrk03)) < fpmin ) {
         //converged!!!
         *density = den;
         *temp = T;
         break;
      }

      preHSE = pre_im1 + .5*(den+den_im1)*drgimh; // HSE pressure
      prim[RHO] = den;
      preEOS = get_pre( prim , &T );
      entEOS = get_entr( prim , &T );
      get_derivs(prim , &T , derivs);
      dpd = derivs[0];
      dpt = derivs[1];
      dsd = derivs[2];
      dst = derivs[3];

      f1 = preEOS - preHSE;
      f2 = entEOS - S_He;
      df11 = dpd - .5*drgimh;
      df12 = dpt;
      df21 = dsd;
      df22 = dst;
      det = df11*df22-df12*df21;
      eoswrk02 = (df22*f1-df12*f2)/det;
      eoswrk03 = (df11*f2-df21*f1)/det;

      dennew = fmin(fmax(.5*den,den - eoswrk02),2.*den);
      tempnew = fmin(fmax(.5*T,T - eoswrk03),2.*T);

      eoswrk01 = fabs((dennew - den)/den) + fabs((tempnew - T)/T);
      den = fmin(1e14,fmax(dennew,1e-11));
      T = fmin(1e13,fmax(tempnew,1e3));
   }
   // if here... not converged :(
   if( i==101 ){ printf("ERROR: Pressure and Entropy failed to balance\n"); }
}
