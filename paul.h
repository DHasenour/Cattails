enum{RHO,PPP,VRR,AAA,XXX};
enum{DDD,TAU,SRR};

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "Reactions/pyna13/network.h"

#define NUM_Q (3+NUM_I+1)
#define NUM_G 2

////////// USEFUL CONSTANTS from pynucastro.constants.constants //////////
// speed of light in vacuum
static const double c_light = 2.99792458e10;  // cm/s
// boltzmann's constant
static const double k_B = 1.3806490000000002e-16;  // erg/K
// planck's constant over 2pi
static const double hbar = 1.0545718176461565e-27;  // erg s
// avogradro's Number
static const double n_A = 6.02214076e+23;  // mol^-1
// radiation constant
static const double a = 7.565733250280009e-15;  //erg/cm^3/K^4
// Stefan-Boltzmann constant
static const double sigma = 5.670374419184432e-05;  //erg/cm^2/K^4/s
// convert eV to erg
static const double ev2erg = 1.6021766339999998e-12;
// convert MeV to eV
static const double MeV2eV = 1.0e6;
// convert MeV to grams
static const double MeV2gr  = (MeV2eV * ev2erg) / (c_light * c_light);
// conversion factor for nuclear energy generation rate
static const double enuc_conv2 = -n_A * c_light * c_light;
// mass of proton
static const double m_p = 1.67262192595e-24;  // g
// mass of neutron
static const double m_n = 1.6749275005600003e-24;  // g
// mass of electron
static const double m_e = 9.1093837139e-28;  // g
// atomic mass unit
static const double m_u = 1.66053906892e-24; // g
// atomic mass unit C18
static const double m_u_C18 = 1.6605390666e-24; // g
// electron charge
static const double q_e = 4.803204712570263e-10;  // erg^1/2 cm^1/2

////////// USEFUL CONSTANTS from helmeos //////////
//stefan-boltzmann constant 
static const double ssol = 5.6704e-5;  //erg cm−2 s−1 K−4
//radiation (density) constant
static const double asol = 4.0 * ssol / c_light;   //erg cm−3 K−4

////////// Other USEFUL CONSTANTS ///////////
//gravitational constant
static const double grav_G = 6.67428e-8;  //cm^3 g-1 s-2 
//solar mass
static const double Msun = 1.9892e33;  //g

struct param_list{

   int Num_R;
   double t_min, t_max;
   double rmin,rmax;
   int NumRepts, NumSnaps, NumChecks;
   int Out_LogTime;

   int LogZoning;
   double LogRadius;
   double LogFractn;
   double LogSmooth;
   int Mesh_Motion, Riemann_Solver;
   double MaxShort, MaxLong;
   int Absorb_BC, Initial_Regrid, rt_flag, rn_flag;

   double ATOL,RTOL;

   int grav_flag, grow_flag;
   double grav_pointmass;
   int grav_e_mode;
   int grav_bal;

   double CFL, PLM;
   double Density_Floor, Pressure_Floor;

   double rt_A,rt_B,rt_C,rt_D;

};

struct domain{

   struct cell * theCells;
   int Nr,Ng;
   double point_mass;

   time_t Wallt_init;
   int rank,size;

   struct param_list theParList;

   double t;
   int count_steps;
   double t_init, t_fin;
   int nrpt;
   int N_rpt;
   int nsnp;
   int N_snp;
   int nchk;
   int N_chk;

   int final_step;

};

struct cell{
   double prim[NUM_Q];
   double cons[NUM_Q];
   double RKcons[NUM_Q];
   double grad[NUM_Q];
   double riph;
   double dr;
   double miph;
   double dm;
   double wiph;
   double T;
   double gradT;
   
};
