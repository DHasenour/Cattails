#ifndef actual_network_H
#define actual_network_H

#include <amrex_bridge.H>

#include <fundamental_constants.H>
#include <network_properties.H>

void actual_network_init();

const std::string network_name = "pynucastro-cxx";

namespace network
{
    extern Array1D<Real, 1, NumSpec> bion;
    extern Array1D<Real, 1, NumSpec> mion;
}

namespace Rates
{

    enum NetworkRates
    {
        k_O16_to_He4_C12 = 1,
        k_Ne20_to_He4_O16 = 2,
        k_Mg24_to_He4_Ne20 = 3,
        k_C12_to_He4_He4_He4 = 4,
        k_He4_C12_to_O16 = 5,
        k_He4_O16_to_Ne20 = 6,
        k_He4_Ne20_to_Mg24 = 7,
        k_C12_C12_to_He4_Ne20 = 8,
        k_C12_O16_to_He4_Mg24 = 9,
        k_O16_O16_to_He4_Si28 = 10,
        k_He4_Ne20_to_C12_C12 = 11,
        k_C12_Ne20_to_He4_Si28 = 12,
        k_He4_Mg24_to_C12_O16 = 13,
        k_He4_Si28_to_C12_Ne20 = 14,
        k_He4_Si28_to_O16_O16 = 15,
        k_He4_He4_He4_to_C12 = 16,
        k_He4_Mg24_to_Si28_removed = 17,
        k_He4_Mg24_to_p_Al27_removed = 18,
        k_p_Al27_to_Si28_removed = 19,
        k_Si28_to_He4_Mg24_removed = 20,
        k_Si28_to_p_Al27_removed = 21,
        k_p_Al27_to_He4_Mg24_removed = 22,
        k_He4_Si28_to_S32_removed = 23,
        k_He4_Si28_to_p_P31_removed = 24,
        k_p_P31_to_S32_removed = 25,
        k_S32_to_He4_Si28_removed = 26,
        k_S32_to_p_P31_removed = 27,
        k_p_P31_to_He4_Si28_removed = 28,
        k_He4_S32_to_Ar36_removed = 29,
        k_He4_S32_to_p_Cl35_removed = 30,
        k_p_Cl35_to_Ar36_removed = 31,
        k_Ar36_to_He4_S32_removed = 32,
        k_Ar36_to_p_Cl35_removed = 33,
        k_p_Cl35_to_He4_S32_removed = 34,
        k_He4_Ar36_to_Ca40_removed = 35,
        k_He4_Ar36_to_p_K39_removed = 36,
        k_p_K39_to_Ca40_removed = 37,
        k_Ca40_to_He4_Ar36_removed = 38,
        k_Ca40_to_p_K39_removed = 39,
        k_p_K39_to_He4_Ar36_removed = 40,
        k_He4_Ca40_to_Ti44_removed = 41,
        k_He4_Ca40_to_p_Sc43_removed = 42,
        k_p_Sc43_to_Ti44_removed = 43,
        k_Ti44_to_He4_Ca40_removed = 44,
        k_Ti44_to_p_Sc43_removed = 45,
        k_p_Sc43_to_He4_Ca40_removed = 46,
        k_He4_Ti44_to_Cr48_removed = 47,
        k_He4_Ti44_to_p_V47_removed = 48,
        k_p_V47_to_Cr48_removed = 49,
        k_Cr48_to_He4_Ti44_removed = 50,
        k_Cr48_to_p_V47_removed = 51,
        k_p_V47_to_He4_Ti44_removed = 52,
        k_He4_Cr48_to_Fe52_removed = 53,
        k_He4_Cr48_to_p_Mn51_removed = 54,
        k_p_Mn51_to_Fe52_removed = 55,
        k_Fe52_to_He4_Cr48_removed = 56,
        k_Fe52_to_p_Mn51_removed = 57,
        k_p_Mn51_to_He4_Cr48_removed = 58,
        k_He4_Fe52_to_Ni56_removed = 59,
        k_He4_Fe52_to_p_Co55_removed = 60,
        k_p_Co55_to_Ni56_removed = 61,
        k_Ni56_to_He4_Fe52_removed = 62,
        k_Ni56_to_p_Co55_removed = 63,
        k_p_Co55_to_He4_Fe52_removed = 64,
        k_Mg24_He4_to_Si28_approx = 65,
        k_Si28_to_Mg24_He4_approx = 66,
        k_Si28_He4_to_S32_approx = 67,
        k_S32_to_Si28_He4_approx = 68,
        k_S32_He4_to_Ar36_approx = 69,
        k_Ar36_to_S32_He4_approx = 70,
        k_Ar36_He4_to_Ca40_approx = 71,
        k_Ca40_to_Ar36_He4_approx = 72,
        k_Ca40_He4_to_Ti44_approx = 73,
        k_Ti44_to_Ca40_He4_approx = 74,
        k_Ti44_He4_to_Cr48_approx = 75,
        k_Cr48_to_Ti44_He4_approx = 76,
        k_Cr48_He4_to_Fe52_approx = 77,
        k_Fe52_to_Cr48_He4_approx = 78,
        k_Fe52_He4_to_Ni56_approx = 79,
        k_Ni56_to_Fe52_He4_approx = 80,
        NumRates = k_Ni56_to_Fe52_He4_approx
    };

    // number of reaclib rates

    const int NrateReaclib = 64;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "O16_to_He4_C12",  // 1,
        "Ne20_to_He4_O16",  // 2,
        "Mg24_to_He4_Ne20",  // 3,
        "C12_to_He4_He4_He4",  // 4,
        "He4_C12_to_O16",  // 5,
        "He4_O16_to_Ne20",  // 6,
        "He4_Ne20_to_Mg24",  // 7,
        "C12_C12_to_He4_Ne20",  // 8,
        "C12_O16_to_He4_Mg24",  // 9,
        "O16_O16_to_He4_Si28",  // 10,
        "He4_Ne20_to_C12_C12",  // 11,
        "C12_Ne20_to_He4_Si28",  // 12,
        "He4_Mg24_to_C12_O16",  // 13,
        "He4_Si28_to_C12_Ne20",  // 14,
        "He4_Si28_to_O16_O16",  // 15,
        "He4_He4_He4_to_C12",  // 16,
        "He4_Mg24_to_Si28_removed",  // 17,
        "He4_Mg24_to_p_Al27_removed",  // 18,
        "p_Al27_to_Si28_removed",  // 19,
        "Si28_to_He4_Mg24_removed",  // 20,
        "Si28_to_p_Al27_removed",  // 21,
        "p_Al27_to_He4_Mg24_removed",  // 22,
        "He4_Si28_to_S32_removed",  // 23,
        "He4_Si28_to_p_P31_removed",  // 24,
        "p_P31_to_S32_removed",  // 25,
        "S32_to_He4_Si28_removed",  // 26,
        "S32_to_p_P31_removed",  // 27,
        "p_P31_to_He4_Si28_removed",  // 28,
        "He4_S32_to_Ar36_removed",  // 29,
        "He4_S32_to_p_Cl35_removed",  // 30,
        "p_Cl35_to_Ar36_removed",  // 31,
        "Ar36_to_He4_S32_removed",  // 32,
        "Ar36_to_p_Cl35_removed",  // 33,
        "p_Cl35_to_He4_S32_removed",  // 34,
        "He4_Ar36_to_Ca40_removed",  // 35,
        "He4_Ar36_to_p_K39_removed",  // 36,
        "p_K39_to_Ca40_removed",  // 37,
        "Ca40_to_He4_Ar36_removed",  // 38,
        "Ca40_to_p_K39_removed",  // 39,
        "p_K39_to_He4_Ar36_removed",  // 40,
        "He4_Ca40_to_Ti44_removed",  // 41,
        "He4_Ca40_to_p_Sc43_removed",  // 42,
        "p_Sc43_to_Ti44_removed",  // 43,
        "Ti44_to_He4_Ca40_removed",  // 44,
        "Ti44_to_p_Sc43_removed",  // 45,
        "p_Sc43_to_He4_Ca40_removed",  // 46,
        "He4_Ti44_to_Cr48_removed",  // 47,
        "He4_Ti44_to_p_V47_removed",  // 48,
        "p_V47_to_Cr48_removed",  // 49,
        "Cr48_to_He4_Ti44_removed",  // 50,
        "Cr48_to_p_V47_removed",  // 51,
        "p_V47_to_He4_Ti44_removed",  // 52,
        "He4_Cr48_to_Fe52_removed",  // 53,
        "He4_Cr48_to_p_Mn51_removed",  // 54,
        "p_Mn51_to_Fe52_removed",  // 55,
        "Fe52_to_He4_Cr48_removed",  // 56,
        "Fe52_to_p_Mn51_removed",  // 57,
        "p_Mn51_to_He4_Cr48_removed",  // 58,
        "He4_Fe52_to_Ni56_removed",  // 59,
        "He4_Fe52_to_p_Co55_removed",  // 60,
        "p_Co55_to_Ni56_removed",  // 61,
        "Ni56_to_He4_Fe52_removed",  // 62,
        "Ni56_to_p_Co55_removed",  // 63,
        "p_Co55_to_He4_Fe52_removed",  // 64,
        "Mg24_He4_to_Si28_approx",  // 65,
        "Si28_to_Mg24_He4_approx",  // 66,
        "Si28_He4_to_S32_approx",  // 67,
        "S32_to_Si28_He4_approx",  // 68,
        "S32_He4_to_Ar36_approx",  // 69,
        "Ar36_to_S32_He4_approx",  // 70,
        "Ar36_He4_to_Ca40_approx",  // 71,
        "Ca40_to_Ar36_He4_approx",  // 72,
        "Ca40_He4_to_Ti44_approx",  // 73,
        "Ti44_to_Ca40_He4_approx",  // 74,
        "Ti44_He4_to_Cr48_approx",  // 75,
        "Cr48_to_Ti44_He4_approx",  // 76,
        "Cr48_He4_to_Fe52_approx",  // 77,
        "Fe52_to_Cr48_He4_approx",  // 78,
        "Fe52_He4_to_Ni56_approx",  // 79,
        "Ni56_to_Fe52_He4_approx"  // 80,
    };

}

#endif
