#include <amrex_bridge.H>
#include <network_properties.H>
#include <actual_rhs.H>

#include <iostream>

extern "C" void derivs(double[],double[],double,double);
void derivs(double y[], double dydx[], double rho, double T) {
  actual_network_init();

  burn_t state;
  state.T = T;
  state.rho = rho;
  for (int n = 0; n < NumSpec; ++n) {state.xn[n] = y[n]*aion[n];}

  Array1D<Real, 1, NumSpec> ydot;
  actual_rhs(state, ydot);

  for (int n = 0; n < NumSpec; ++n) {dydx[n] = ydot(n+1);}
}

extern "C" void jacobn(double[],double[][NumSpec],double,double);
void jacobn(double y[], double dfdy[][NumSpec], double rho, double T) {
  actual_network_init();

  burn_t state;
  state.T = T;
  state.rho = rho;
  for (int n = 0; n < NumSpec; ++n) {state.xn[n] = y[n]*aion[n];}

  MathArray2D<1, NumSpec, 1, NumSpec> jac;
  actual_jac(state, jac);

  for (int n = 0; n < NumSpec; ++n) {
    for (int m = 0; m < NumSpec; ++m) {
      dfdy[n][m] = jac(n+1,m+1);
    }
  }
}
