
#include "../../paul.h"

void getcons_( double * , double * , double * , double * , double * , double * , double * );
void getprim_( double * , double * , double * , double * , double * , double * );
void getambient_( double * , double * , double * , double * , double * , double * , double * , double * , double * );
double get_dV( double , double );
void calc_prim( struct domain * );
void boundary( struct domain * );
void exchangeData( struct domain * );

double RTOL = 1.e-7;
double ATOL = 1.e-7;

void derivs(double[],double[],double,double);
void jacobn(double[],double[][NUM_I],double,double);

void ludcmp(double a[NUM_I][NUM_I], int indx[], double *d) {
  double TINY = 1.0e-20;
  int i,imax=0,j,k;
  double big,dum,sum,temp;
  double vv[NUM_I];

  *d=1.0;
  for (i=0;i<NUM_I;i++) {
    big=0.0;
    for (j=0;j<NUM_I;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    //if (big == 0.0) printf("Singular matrix in routine ludcmp\n");
    vv[i]=1.0/big;
  }
  for (j=0;j<NUM_I;j++) {
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<NUM_I;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<NUM_I;k++) {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != NUM_I-1) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<NUM_I;i++) a[i][j] *= dum;
    }
  }
}
void lubksb(double a[NUM_I][NUM_I], int indx[], double b[]) {
  int i,ii=-1,ip,j;
  double sum;

  for (i=0;i<NUM_I;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii>=0)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=NUM_I-1;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<NUM_I;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

void odeint( double ystart[] , double dt , double rho , double T ){
  double d , x;
  double dfdy[NUM_I][NUM_I] , a[NUM_I][NUM_I];
  double dydx[NUM_I] , y1[NUM_I] , y2[NUM_I] , y3[NUM_I] , y4[NUM_I] , yout[NUM_I] , yscal[NUM_I];
  int i , j , indx[NUM_I];
  int nok , nbad , norder , stepMAX , stepNOW;

  nok = 0;
  nbad = 0;
  norder = 2;
  stepMAX = 10000000;
  stepNOW = 0;

  /////// Initialize hydrostatic specific total energy: e_cons = e_bind + e_int ////////
  double e_cons , e_int , mass_frac[NUM_I] , pre , s , gam , cs , cv;
  for( int i=0;i<NUM_I;i++ ){ mass_frac[i]=ystart[i]*aion[i]; }
  getambient_( &rho , &T , mass_frac , &pre , &s , &gam , &e_int , &cs , &cv );
  e_cons = e_int;
  for( int i=0;i<NUM_I;i++ ){ e_cons+=mass_frac[i]*EBIND[i]; }

  double hscale;
  double del_x = dt;
  x = 0.0;
  while ( x<dt ){
    // limit max steps to not forever run 
    stepNOW += 1;
    if (stepNOW>stepMAX) {
        printf("ERROR: Max steps reached in ODEINT\n");
        exit(0);
    }

    // Initialize y1 and y2
    for (i=0;i<NUM_I;i++) {
      y1[i]=ystart[i];
      y2[i]=ystart[i];
      y3[i]=ystart[i];
      y4[i]=ystart[i];
    }
    // Desired accuracy for step. Mostly Constant fractional errors
    for (i=0;i<NUM_I;i++) {yscal[i]=(RTOL*fabs(ystart[i])+ATOL);}

    // Temperature Update //
    e_int = e_cons;
    for( i=0;i<NUM_I;i++ ){
      mass_frac[i]=ystart[i]*aion[i];
      e_int-=mass_frac[i]*EBIND[i];
    }
    getprim_( &rho , &e_int , mass_frac , &pre , &T , &cs );

//////////////////////////////////////////////////////////////////////////////////
    if (norder==2) {
      /////////// Full Step /////////////
      jacobn(y1,dfdy,rho,T);
      derivs(y1,dydx,rho,T);

      for (i=0;i<NUM_I;i++) {
        for (j=0;j<NUM_I;j++) {a[i][j] = -del_x*dfdy[i][j];}
        a[i][i] += 1.0;
      }
      ludcmp(a,indx,&d);
      for (i=0;i<NUM_I;i++) {yout[i]=del_x*dydx[i];}
      lubksb(a,indx,yout);
      for (i=0;i<NUM_I;i++) {y1[i]+=yout[i];}

      for (i=0;i<NUM_I;i++) {
        for (j=0;j<NUM_I;j++) {a[i][j] = -.5*del_x*dfdy[i][j];}
        a[i][i] += 1.0;
      }
      ludcmp(a,indx,&d);
      for (i=0;i<NUM_I;i++) {yout[i]=.5*del_x*dydx[i];}
      lubksb(a,indx,yout);
      for (i=0;i<NUM_I;i++) {y2[i]+=yout[i];}
      
      ///////// Half Steps ///////////////
      jacobn(y2,dfdy,rho,T);
      derivs(y2,dydx,rho,T);
      
      for (i=0;i<NUM_I;i++) {
        for (j=0;j<NUM_I;j++) {a[i][j] = -.5*del_x*dfdy[i][j];}
        a[i][i] += 1.0;
      }
      ludcmp(a,indx,&d);
      for (i=0;i<NUM_I;i++) {yout[i]=.5*del_x*dydx[i];}
      lubksb(a,indx,yout);
      for (i=0;i<NUM_I;i++) {y2[i]+=yout[i];}

      ////////// Adaptive Step ////////////
      for (i=0;i<NUM_I;i++) {yout[i] = fabs(yscal[i]/(y2[i]-y1[i]));}
      hscale = yout[1];
      for (i=1;i<NUM_I;i++) {
          if (yout[i]<hscale) {hscale=yout[i];}
      }
      hscale = sqrt(hscale);
      if(hscale != hscale){hscale = 0.5;}  //assume a nan is a failed state

      if (hscale >= 1.) {
        //succesful step
        for (i=0;i<NUM_I;i++) {ystart[i]=2.*y2[i]-y1[i];}
        x += del_x;
        del_x = del_x*hscale;
        if(x+del_x>=dt) {del_x = dt-x;}  //Special last step
        nok+=1;
      }
      else {
        //failed step 
        del_x = del_x*hscale*0.5; //0.5 is an arb choice
        if(x+del_x>=dt) {del_x = dt-x;}  //Special last step
        nbad+=1;
        norder=4;
      }
    }

    else {
      /////////// Full Step /////////////
      jacobn(y1,dfdy,rho,T);
      derivs(y1,dydx,rho,T);

      for (i=0;i<NUM_I;i++) {
        for (j=0;j<NUM_I;j++) {a[i][j] = -del_x*dfdy[i][j];}
        a[i][i] += 1.0;
      }
      ludcmp(a,indx,&d);
      for (i=0;i<NUM_I;i++) {yout[i]=del_x*dydx[i];}
      lubksb(a,indx,yout);
      for (i=0;i<NUM_I;i++) {y1[i]+=yout[i];}

      for (i=0;i<NUM_I;i++) {
        for (j=0;j<NUM_I;j++) {a[i][j] = -.5*del_x*dfdy[i][j];}
        a[i][i] += 1.0;
      }
      ludcmp(a,indx,&d);
      for (i=0;i<NUM_I;i++) {yout[i]=.5*del_x*dydx[i];}
      lubksb(a,indx,yout);
      for (i=0;i<NUM_I;i++) {y2[i]+=yout[i];}

      for (i=0;i<NUM_I;i++) {
        for (j=0;j<NUM_I;j++) {a[i][j] = -(1./3.)*del_x*dfdy[i][j];}
        a[i][i] += 1.0;
      }
      ludcmp(a,indx,&d);
      for (i=0;i<NUM_I;i++) {yout[i]=(1./3.)*del_x*dydx[i];}
      lubksb(a,indx,yout);
      for (i=0;i<NUM_I;i++) {y3[i]+=yout[i];}

      for (i=0;i<NUM_I;i++) {
        for (j=0;j<NUM_I;j++) {a[i][j] = -.25*del_x*dfdy[i][j];}
        a[i][i] += 1.0;
      }
      ludcmp(a,indx,&d);
      for (i=0;i<NUM_I;i++) {yout[i]=.25*del_x*dydx[i];}
      lubksb(a,indx,yout);
      for (i=0;i<NUM_I;i++) {y4[i]+=yout[i];}
      
      ///////// Half Steps ///////////////
      jacobn(y2,dfdy,rho,T);
      derivs(y2,dydx,rho,T);
      
      for (i=0;i<NUM_I;i++) {
        for (j=0;j<NUM_I;j++) {a[i][j] = -.5*del_x*dfdy[i][j];}
        a[i][i] += 1.0;
      }
      ludcmp(a,indx,&d);
      for (i=0;i<NUM_I;i++) {yout[i]=.5*del_x*dydx[i];}
      lubksb(a,indx,yout);
      for (i=0;i<NUM_I;i++) {y2[i]+=yout[i];}

      ///////// THIRD STEPS //////////////
      jacobn(y3,dfdy,rho,T);
      derivs(y3,dydx,rho,T);
      
      for (i=0;i<NUM_I;i++) {
        for (j=0;j<NUM_I;j++) {a[i][j] = -(1./3.)*del_x*dfdy[i][j];}
        a[i][i] += 1.0;
      }
      ludcmp(a,indx,&d);
      for (i=0;i<NUM_I;i++) {yout[i]=(1./3.)*del_x*dydx[i];}
      lubksb(a,indx,yout);
      for (i=0;i<NUM_I;i++) {y3[i]+=yout[i];}

      jacobn(y3,dfdy,rho,T);
      derivs(y3,dydx,rho,T);
      
      for (i=0;i<NUM_I;i++) {
        for (j=0;j<NUM_I;j++) {a[i][j] = -(1./3.)*del_x*dfdy[i][j];}
        a[i][i] += 1.0;
      }
      ludcmp(a,indx,&d);
      for (i=0;i<NUM_I;i++) {yout[i]=(1./3.)*del_x*dydx[i];}
      lubksb(a,indx,yout);
      for (i=0;i<NUM_I;i++) {y3[i]+=yout[i];}

      ///////// Fourth Steps ///////////////
      jacobn(y4,dfdy,rho,T);
      derivs(y4,dydx,rho,T);
      
      for (i=0;i<NUM_I;i++) {
        for (j=0;j<NUM_I;j++) {a[i][j] = -.25*del_x*dfdy[i][j];}
        a[i][i] += 1.0;
      }
      ludcmp(a,indx,&d);
      for (i=0;i<NUM_I;i++) {yout[i]=.25*del_x*dydx[i];}
      lubksb(a,indx,yout);
      for (i=0;i<NUM_I;i++) {y4[i]+=yout[i];}

      jacobn(y4,dfdy,rho,T);
      derivs(y4,dydx,rho,T);
      
      for (i=0;i<NUM_I;i++) {
        for (j=0;j<NUM_I;j++) {a[i][j] = -.25*del_x*dfdy[i][j];}
        a[i][i] += 1.0;
      }
      ludcmp(a,indx,&d);
      for (i=0;i<NUM_I;i++) {yout[i]=.25*del_x*dydx[i];}
      lubksb(a,indx,yout);
      for (i=0;i<NUM_I;i++) {y4[i]+=yout[i];}

      jacobn(y4,dfdy,rho,T);
      derivs(y4,dydx,rho,T);
      
      for (i=0;i<NUM_I;i++) {
        for (j=0;j<NUM_I;j++) {a[i][j] = -.25*del_x*dfdy[i][j];}
        a[i][i] += 1.0;
      }
      ludcmp(a,indx,&d);
      for (i=0;i<NUM_I;i++) {yout[i]=.25*del_x*dydx[i];}
      lubksb(a,indx,yout);
      for (i=0;i<NUM_I;i++) {y4[i]+=yout[i];}

      ////////// Adaptive Step ////////////
      for (i=0;i<NUM_I;i++) {yout[i] = fabs(yscal[i]/((8./3.)*y4[i]-4.5*y3[i]+2.*y2[i]-(1./6.)*y1[i]));}
      hscale = yout[1];
      for (i=1;i<NUM_I;i++) {
          if (yout[i]<hscale) {hscale=yout[i];}
      }
      hscale = pow(hscale,.25);
      if(hscale != hscale){hscale = 0.5;}  //assume a nan is a failed state

      if (hscale >= 1.) {
        //succesful step
        for (i=0;i<NUM_I;i++) {ystart[i]=((32./3.)*y4[i]-13.5*y3[i]+4.*y2[i]-(1./6.)*y1[i]);}
        x += del_x;
        del_x = del_x*hscale;
        if(x+del_x>=dt) {del_x = dt-x;}  //Special last step
        nok+=1;
      }
      else {
        //failed step 
        del_x = del_x*hscale*0.5; //0.5 is an arb choice
        if(x+del_x>=dt) {del_x = dt-x;}  //Special last step
        nbad+=1;
      }
    }
  }
  //printf("nok = %d  nbad = %d\n",nok,nbad);
}

void reacstep( struct domain * theDomain , double dt ) {
  double rho , T , ystart[NUM_I];
  double Pp , etot , stot , cs;
  struct cell * theCells = theDomain->theCells;
  int Nr = theDomain->Nr;
  int i, j;
  double rm , rp , dV; 

  /*
  double array[N];
  double ptot,stot,gam,etot,cs,cv;
  double global_tsound = 1e199, local_tsound;

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Burning limiter: min sound-crossing time //
      for(int i=2;i<N+2;i++){
          for (int j=0;j<NUM_I;j++) {ystart[j]=cells[i].prim[3+j];} // really is mass fraction x[j]
          rho = cells[i].prim[DEN];
          T = cells[i].T;
          getambient_(&rho,&T,ystart,&ptot,&stot,&gam,&etot,&cs,&cv);
          array[(i-2)] = cells[i].dV/cs;
      }
      local_tsound = MIN(array,N);
      MPI_Allreduce(&local_tsound,&global_tsound,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  */

  //Loop over cells
  for( i=0 ; i<Nr ; ++i ) {
    struct cell * c = theCells+i;
    for( j=0 ; j<NUM_I ; ++j ){ ystart[j] = c->prim[3+j]; }
    rho = c->prim[RHO];
    Pp = c->prim[PPP];
    T = c->T; //a guess
    getcons_( &rho , &Pp , ystart , &etot , &stot , &T , &cs );
    for( j=0 ; j<NUM_I ; ++j ){ ystart[j] = ystart[j]/aion[j]; }

    //if( T>1.1e7 & c->riph>1e7 ) { odeint( ystart , dt , rho , T ); }
    if( T>1.1e7 ) { odeint( ystart , dt , rho , T ); }
    //odeint( ystart , dt , rho , T );
      
    rp = c->riph;
    rm = rp-c->dr;
    dV = get_dV( rp , rm );
    for( j=0 ; j<NUM_I ; ++j ){ c->cons[3+j] = ystart[j]*rho*dV*aion[j];}
  }
  calc_prim( theDomain );
  boundary( theDomain );
  exchangeData( theDomain );
}
