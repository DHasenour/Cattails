#ifndef actual_rhs_H
#define actual_rhs_H

#include <amrex_bridge.H>

#include <actual_network.H>
#include <burn_type.H>
#include <screen.H>
#include <reaclib_rates.H>

using namespace Species;
using namespace Rates;


template<class T>
inline
void ener_gener_rate(T const& dydt, Real& enuc)
{

    // Computes the instantaneous energy generation rate (from the nuclei)

    // This is basically e = m c**2

    enuc = 0.0_rt;

    for (int n = 1; n <= NumSpec; ++n) {
        enuc += dydt(n) * network::mion(n);
    }

    enuc *= C::enuc_conv2;
}


template <int do_T_derivatives, typename T>
inline
void evaluate_rates(const burn_t& state, T& rate_eval) {

    // create molar fractions

    Array1D<Real, 1, NumSpec> Y;
    for (int n = 1; n <= NumSpec; ++n) {
        Y(n) = state.xn[n-1] * aion_inv[n-1];
    }

    // Calculate Reaclib rates

    plasma_state_t pstate{};
    fill_plasma_state(pstate, state.T, state.rho, Y);

    tf_t tfactors = evaluate_tfactors(state.T);

    fill_reaclib_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // Evaluating screening factors

    Real scor;
    [[maybe_unused]] Real scor2;


    // Fill approximate rates

    fill_approx_rates<do_T_derivatives, T>(tfactors, state.rho, Y, rate_eval);

}

inline
void rhs_nuc(const burn_t& state,
             Array1D<Real, 1, NumSpec>& ydot_nuc,
             const Array1D<Real, 1, NumSpec>& Y,
             const Array1D<Real, 1, NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(He4) =
        (-screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12)*Y(O16)) +
        (-screened_rates(k_He4_O16_to_Ne20)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16)*Y(Ne20)) +
        (-screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20)*Y(Mg24)) +
        (0.5*screened_rates(k_C12_C12_to_He4_Ne20)*std::pow(Y(C12), 2)*state.rho + -screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*Y(Ne20)*state.rho) +
        (screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*Y(Mg24)*state.rho) +
        (0.5*screened_rates(k_O16_O16_to_He4_Si28)*std::pow(Y(O16), 2)*state.rho + -screened_rates(k_He4_Si28_to_O16_O16)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*Y(Ne20)*state.rho + -screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*Y(Si28)*state.rho) +
        (-0.5*screened_rates(k_He4_He4_He4_to_C12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2) + 3.0*screened_rates(k_C12_to_He4_He4_He4)*Y(C12)) +
        (-screened_rates(k_Mg24_He4_to_Si28_approx)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_Si28_to_Mg24_He4_approx)*Y(Si28)) +
        (-screened_rates(k_Si28_He4_to_S32_approx)*Y(He4)*Y(Si28)*state.rho + screened_rates(k_S32_to_Si28_He4_approx)*Y(S32)) +
        (-screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*Y(S32)*state.rho + screened_rates(k_Ar36_to_S32_He4_approx)*Y(Ar36)) +
        (-screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*Y(He4)*state.rho + screened_rates(k_Ca40_to_Ar36_He4_approx)*Y(Ca40)) +
        (-screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*Y(He4)*state.rho + screened_rates(k_Ti44_to_Ca40_He4_approx)*Y(Ti44)) +
        (-screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*Y(Ti44)*state.rho + screened_rates(k_Cr48_to_Ti44_He4_approx)*Y(Cr48)) +
        (-screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_Fe52_to_Cr48_He4_approx)*Y(Fe52)) +
        (-screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Ni56_to_Fe52_He4_approx)*Y(Ni56));

    ydot_nuc(C12) =
        (-screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12)*Y(O16)) +
        (-screened_rates(k_C12_C12_to_He4_Ne20)*std::pow(Y(C12), 2)*state.rho + 2.0*screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho + screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*Y(Ne20)*state.rho + screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*Y(Si28)*state.rho) +
        (0.16666666666666667*screened_rates(k_He4_He4_He4_to_C12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2) + -screened_rates(k_C12_to_He4_He4_He4)*Y(C12));

    ydot_nuc(O16) =
        (screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho + -screened_rates(k_O16_to_He4_C12)*Y(O16)) +
        (-screened_rates(k_He4_O16_to_Ne20)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16)*Y(Ne20)) +
        (-screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho + screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_O16_O16_to_He4_Si28)*std::pow(Y(O16), 2)*state.rho + 2.0*screened_rates(k_He4_Si28_to_O16_O16)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(Ne20) =
        (screened_rates(k_He4_O16_to_Ne20)*Y(He4)*Y(O16)*state.rho + -screened_rates(k_Ne20_to_He4_O16)*Y(Ne20)) +
        (-screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20)*Y(Mg24)) +
        (0.5*screened_rates(k_C12_C12_to_He4_Ne20)*std::pow(Y(C12), 2)*state.rho + -screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*Y(Ne20)*state.rho + screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(Mg24) =
        (screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*Y(Ne20)*state.rho + -screened_rates(k_Mg24_to_He4_Ne20)*Y(Mg24)) +
        (screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_Mg24_He4_to_Si28_approx)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_Si28_to_Mg24_He4_approx)*Y(Si28));

    ydot_nuc(Si28) =
        (0.5*screened_rates(k_O16_O16_to_He4_Si28)*std::pow(Y(O16), 2)*state.rho + -screened_rates(k_He4_Si28_to_O16_O16)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*Y(Ne20)*state.rho + -screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_Mg24_He4_to_Si28_approx)*Y(He4)*Y(Mg24)*state.rho + -screened_rates(k_Si28_to_Mg24_He4_approx)*Y(Si28)) +
        (-screened_rates(k_Si28_He4_to_S32_approx)*Y(He4)*Y(Si28)*state.rho + screened_rates(k_S32_to_Si28_He4_approx)*Y(S32));

    ydot_nuc(S32) =
        (screened_rates(k_Si28_He4_to_S32_approx)*Y(He4)*Y(Si28)*state.rho + -screened_rates(k_S32_to_Si28_He4_approx)*Y(S32)) +
        (-screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*Y(S32)*state.rho + screened_rates(k_Ar36_to_S32_He4_approx)*Y(Ar36));

    ydot_nuc(Ar36) =
        (screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*Y(S32)*state.rho + -screened_rates(k_Ar36_to_S32_He4_approx)*Y(Ar36)) +
        (-screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*Y(He4)*state.rho + screened_rates(k_Ca40_to_Ar36_He4_approx)*Y(Ca40));

    ydot_nuc(Ca40) =
        (screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*Y(He4)*state.rho + -screened_rates(k_Ca40_to_Ar36_He4_approx)*Y(Ca40)) +
        (-screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*Y(He4)*state.rho + screened_rates(k_Ti44_to_Ca40_He4_approx)*Y(Ti44));

    ydot_nuc(Ti44) =
        (screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*Y(He4)*state.rho + -screened_rates(k_Ti44_to_Ca40_He4_approx)*Y(Ti44)) +
        (-screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*Y(Ti44)*state.rho + screened_rates(k_Cr48_to_Ti44_He4_approx)*Y(Cr48));

    ydot_nuc(Cr48) =
        (screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*Y(Ti44)*state.rho + -screened_rates(k_Cr48_to_Ti44_He4_approx)*Y(Cr48)) +
        (-screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_Fe52_to_Cr48_He4_approx)*Y(Fe52));

    ydot_nuc(Fe52) =
        (screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_Fe52_to_Cr48_He4_approx)*Y(Fe52)) +
        (-screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Ni56_to_Fe52_He4_approx)*Y(Ni56));

    ydot_nuc(Ni56) =
        (screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_Ni56_to_Fe52_He4_approx)*Y(Ni56));

}


inline
void actual_rhs (burn_t& state, Array1D<Real, 1, NumSpec>& ydot)
{
    for (int i = 1; i <= NumSpec; ++i) {
        ydot(i) = 0.0_rt;
    }


    // Set molar abundances
    Array1D<Real, 1, NumSpec> Y;
    for (int i = 1; i <= NumSpec; ++i) {
        Y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // build the rates

    rate_t rate_eval;

    constexpr int do_T_derivatives = 0;
    evaluate_rates<do_T_derivatives, rate_t>(state, rate_eval);

    rhs_nuc(state, ydot, Y, rate_eval.screened_rates);

}


template<class MatrixType>
inline
void jac_nuc(const burn_t& state,
             MatrixType& jac,
             const Array1D<Real, 1, NumSpec>& Y,
             const Array1D<Real, 1, NumRates>& screened_rates)
{

    Real scratch;

    scratch = -screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*state.rho - screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*state.rho - screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(Cr48)*state.rho - screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(Fe52)*state.rho - screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho - 1.5*screened_rates(k_He4_He4_He4_to_C12)*std::pow(Y(He4), 2)*std::pow(state.rho, 2) - screened_rates(k_He4_Mg24_to_C12_O16)*Y(Mg24)*state.rho - screened_rates(k_He4_Ne20_to_C12_C12)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_Mg24)*Y(Ne20)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(O16)*state.rho - screened_rates(k_He4_Si28_to_C12_Ne20)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_O16_O16)*Y(Si28)*state.rho - screened_rates(k_Mg24_He4_to_Si28_approx)*Y(Mg24)*state.rho - screened_rates(k_S32_He4_to_Ar36_approx)*Y(S32)*state.rho - screened_rates(k_Si28_He4_to_S32_approx)*Y(Si28)*state.rho - screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(Ti44)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20)*Y(C12)*state.rho + screened_rates(k_C12_Ne20_to_He4_Si28)*Y(Ne20)*state.rho + screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho + 3.0*screened_rates(k_C12_to_He4_He4_He4) - screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho;
    jac.set(He4, C12, scratch);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(He4)*state.rho + 1.0*screened_rates(k_O16_O16_to_He4_Si28)*Y(O16)*state.rho + screened_rates(k_O16_to_He4_C12);
    jac.set(He4, O16, scratch);

    scratch = screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*state.rho - screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*state.rho + screened_rates(k_Ne20_to_He4_O16);
    jac.set(He4, Ne20, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*state.rho - screened_rates(k_Mg24_He4_to_Si28_approx)*Y(He4)*state.rho + screened_rates(k_Mg24_to_He4_Ne20);
    jac.set(He4, Mg24, scratch);

    scratch = -screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_O16_O16)*Y(He4)*state.rho - screened_rates(k_Si28_He4_to_S32_approx)*Y(He4)*state.rho + screened_rates(k_Si28_to_Mg24_He4_approx);
    jac.set(He4, Si28, scratch);

    scratch = -screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*state.rho + screened_rates(k_S32_to_Si28_He4_approx);
    jac.set(He4, S32, scratch);

    scratch = -screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(He4)*state.rho + screened_rates(k_Ar36_to_S32_He4_approx);
    jac.set(He4, Ar36, scratch);

    scratch = -screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(He4)*state.rho + screened_rates(k_Ca40_to_Ar36_He4_approx);
    jac.set(He4, Ca40, scratch);

    scratch = -screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*state.rho + screened_rates(k_Ti44_to_Ca40_He4_approx);
    jac.set(He4, Ti44, scratch);

    scratch = -screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(He4)*state.rho + screened_rates(k_Cr48_to_Ti44_He4_approx);
    jac.set(He4, Cr48, scratch);

    scratch = -screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(He4)*state.rho + screened_rates(k_Fe52_to_Cr48_He4_approx);
    jac.set(He4, Fe52, scratch);

    scratch = screened_rates(k_Ni56_to_Fe52_He4_approx);
    jac.set(He4, Ni56, scratch);

    scratch = -screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho + 0.5*screened_rates(k_He4_He4_He4_to_C12)*std::pow(Y(He4), 2)*std::pow(state.rho, 2) + screened_rates(k_He4_Mg24_to_C12_O16)*Y(Mg24)*state.rho + 2.0*screened_rates(k_He4_Ne20_to_C12_C12)*Y(Ne20)*state.rho + screened_rates(k_He4_Si28_to_C12_Ne20)*Y(Si28)*state.rho;
    jac.set(C12, He4, scratch);

    scratch = -2.0*screened_rates(k_C12_C12_to_He4_Ne20)*Y(C12)*state.rho - screened_rates(k_C12_Ne20_to_He4_Si28)*Y(Ne20)*state.rho - screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho - screened_rates(k_C12_to_He4_He4_He4) - screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho;
    jac.set(C12, C12, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho + screened_rates(k_O16_to_He4_C12);
    jac.set(C12, O16, scratch);

    scratch = -screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*state.rho + 2.0*screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*state.rho;
    jac.set(C12, Ne20, scratch);

    scratch = screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*state.rho;
    jac.set(C12, Mg24, scratch);

    scratch = screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*state.rho;
    jac.set(C12, Si28, scratch);

    scratch = screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho + screened_rates(k_He4_Mg24_to_C12_O16)*Y(Mg24)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(O16)*state.rho + 2.0*screened_rates(k_He4_Si28_to_O16_O16)*Y(Si28)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho + screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(He4)*state.rho - 2.0*screened_rates(k_O16_O16_to_He4_Si28)*Y(O16)*state.rho - screened_rates(k_O16_to_He4_C12);
    jac.set(O16, O16, scratch);

    scratch = screened_rates(k_Ne20_to_He4_O16);
    jac.set(O16, Ne20, scratch);

    scratch = screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*state.rho;
    jac.set(O16, Mg24, scratch);

    scratch = 2.0*screened_rates(k_He4_Si28_to_O16_O16)*Y(He4)*state.rho;
    jac.set(O16, Si28, scratch);

    scratch = -screened_rates(k_He4_Ne20_to_C12_C12)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_Mg24)*Y(Ne20)*state.rho + screened_rates(k_He4_O16_to_Ne20)*Y(O16)*state.rho + screened_rates(k_He4_Si28_to_C12_Ne20)*Y(Si28)*state.rho;
    jac.set(Ne20, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20)*Y(C12)*state.rho - screened_rates(k_C12_Ne20_to_He4_Si28)*Y(Ne20)*state.rho;
    jac.set(Ne20, C12, scratch);

    scratch = screened_rates(k_He4_O16_to_Ne20)*Y(He4)*state.rho;
    jac.set(Ne20, O16, scratch);

    scratch = -screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*state.rho - screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*state.rho - screened_rates(k_Ne20_to_He4_O16);
    jac.set(Ne20, Ne20, scratch);

    scratch = screened_rates(k_Mg24_to_He4_Ne20);
    jac.set(Ne20, Mg24, scratch);

    scratch = screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*state.rho;
    jac.set(Ne20, Si28, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_C12_O16)*Y(Mg24)*state.rho + screened_rates(k_He4_Ne20_to_Mg24)*Y(Ne20)*state.rho - screened_rates(k_Mg24_He4_to_Si28_approx)*Y(Mg24)*state.rho;
    jac.set(Mg24, He4, scratch);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho;
    jac.set(Mg24, C12, scratch);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho;
    jac.set(Mg24, O16, scratch);

    scratch = screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*state.rho;
    jac.set(Mg24, Ne20, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*state.rho - screened_rates(k_Mg24_He4_to_Si28_approx)*Y(He4)*state.rho - screened_rates(k_Mg24_to_He4_Ne20);
    jac.set(Mg24, Mg24, scratch);

    scratch = screened_rates(k_Si28_to_Mg24_He4_approx);
    jac.set(Mg24, Si28, scratch);

    scratch = -screened_rates(k_He4_Si28_to_C12_Ne20)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_O16_O16)*Y(Si28)*state.rho + screened_rates(k_Mg24_He4_to_Si28_approx)*Y(Mg24)*state.rho - screened_rates(k_Si28_He4_to_S32_approx)*Y(Si28)*state.rho;
    jac.set(Si28, He4, scratch);

    scratch = screened_rates(k_C12_Ne20_to_He4_Si28)*Y(Ne20)*state.rho;
    jac.set(Si28, C12, scratch);

    scratch = 1.0*screened_rates(k_O16_O16_to_He4_Si28)*Y(O16)*state.rho;
    jac.set(Si28, O16, scratch);

    scratch = screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*state.rho;
    jac.set(Si28, Ne20, scratch);

    scratch = screened_rates(k_Mg24_He4_to_Si28_approx)*Y(He4)*state.rho;
    jac.set(Si28, Mg24, scratch);

    scratch = -screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_O16_O16)*Y(He4)*state.rho - screened_rates(k_Si28_He4_to_S32_approx)*Y(He4)*state.rho - screened_rates(k_Si28_to_Mg24_He4_approx);
    jac.set(Si28, Si28, scratch);

    scratch = screened_rates(k_S32_to_Si28_He4_approx);
    jac.set(Si28, S32, scratch);

    scratch = -screened_rates(k_S32_He4_to_Ar36_approx)*Y(S32)*state.rho + screened_rates(k_Si28_He4_to_S32_approx)*Y(Si28)*state.rho;
    jac.set(S32, He4, scratch);

    scratch = screened_rates(k_Si28_He4_to_S32_approx)*Y(He4)*state.rho;
    jac.set(S32, Si28, scratch);

    scratch = -screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*state.rho - screened_rates(k_S32_to_Si28_He4_approx);
    jac.set(S32, S32, scratch);

    scratch = screened_rates(k_Ar36_to_S32_He4_approx);
    jac.set(S32, Ar36, scratch);

    scratch = -screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*state.rho + screened_rates(k_S32_He4_to_Ar36_approx)*Y(S32)*state.rho;
    jac.set(Ar36, He4, scratch);

    scratch = screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*state.rho;
    jac.set(Ar36, S32, scratch);

    scratch = -screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(He4)*state.rho - screened_rates(k_Ar36_to_S32_He4_approx);
    jac.set(Ar36, Ar36, scratch);

    scratch = screened_rates(k_Ca40_to_Ar36_He4_approx);
    jac.set(Ar36, Ca40, scratch);

    scratch = screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*state.rho - screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*state.rho;
    jac.set(Ca40, He4, scratch);

    scratch = screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(He4)*state.rho;
    jac.set(Ca40, Ar36, scratch);

    scratch = -screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(He4)*state.rho - screened_rates(k_Ca40_to_Ar36_He4_approx);
    jac.set(Ca40, Ca40, scratch);

    scratch = screened_rates(k_Ti44_to_Ca40_He4_approx);
    jac.set(Ca40, Ti44, scratch);

    scratch = screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*state.rho - screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(Ti44)*state.rho;
    jac.set(Ti44, He4, scratch);

    scratch = screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(He4)*state.rho;
    jac.set(Ti44, Ca40, scratch);

    scratch = -screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*state.rho - screened_rates(k_Ti44_to_Ca40_He4_approx);
    jac.set(Ti44, Ti44, scratch);

    scratch = screened_rates(k_Cr48_to_Ti44_He4_approx);
    jac.set(Ti44, Cr48, scratch);

    scratch = -screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(Cr48)*state.rho + screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(Ti44)*state.rho;
    jac.set(Cr48, He4, scratch);

    scratch = screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*state.rho;
    jac.set(Cr48, Ti44, scratch);

    scratch = -screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(He4)*state.rho - screened_rates(k_Cr48_to_Ti44_He4_approx);
    jac.set(Cr48, Cr48, scratch);

    scratch = screened_rates(k_Fe52_to_Cr48_He4_approx);
    jac.set(Cr48, Fe52, scratch);

    scratch = screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(Cr48)*state.rho - screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(Fe52)*state.rho;
    jac.set(Fe52, He4, scratch);

    scratch = screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(He4)*state.rho;
    jac.set(Fe52, Cr48, scratch);

    scratch = -screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(He4)*state.rho - screened_rates(k_Fe52_to_Cr48_He4_approx);
    jac.set(Fe52, Fe52, scratch);

    scratch = screened_rates(k_Ni56_to_Fe52_He4_approx);
    jac.set(Fe52, Ni56, scratch);

    scratch = screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(Fe52)*state.rho;
    jac.set(Ni56, He4, scratch);

    scratch = screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(He4)*state.rho;
    jac.set(Ni56, Fe52, scratch);

    scratch = -screened_rates(k_Ni56_to_Fe52_He4_approx);
    jac.set(Ni56, Ni56, scratch);


}



template<class MatrixType>
inline
void actual_jac(const burn_t& state, MatrixType& jac)
{

    // Set molar abundances
    Array1D<Real, 1, NumSpec> Y;
    for (int i = 1; i <= NumSpec; ++i) {
        Y(i) = state.xn[i-1] * aion_inv[i-1];
    }


    jac.zero();

    rate_derivs_t rate_eval;

    constexpr int do_T_derivatives = 0;
    evaluate_rates<do_T_derivatives, rate_derivs_t>(state, rate_eval);

    // Species Jacobian elements with respect to other species

    jac_nuc(state, jac, Y, rate_eval.screened_rates);

}

#endif
