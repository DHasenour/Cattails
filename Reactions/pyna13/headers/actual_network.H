#ifndef actual_network_H
#define actual_network_H

#include <amrex_bridge.H>

#include <fundamental_constants.H>
#include <network_properties.H>

inline void actual_network_init() {}

const std::string network_name = "pynucastro-cxx";

namespace network
{

    // nuclei masses (in grams)

    inline Array1D<Real, 1, NumSpec> mion{
        6.646479071584587e-24_rt,  // He4
        1.99264687992e-23_rt,  // C12
        2.6560180592333686e-23_rt,  // O16
        3.3198227947612416e-23_rt,  // Ne20
        3.9828098739467446e-23_rt,  // Mg24
        4.6456779473820677e-23_rt,  // Si28
        5.309087322384128e-23_rt,  // S32
        5.972551377884467e-23_rt,  // Ar36
        6.635944331004904e-23_rt,  // Ca40
        7.299678247096977e-23_rt,  // Ti44
        7.962953983065421e-23_rt,  // Cr48
        8.626187166893794e-23_rt,  // Fe52
        9.289408870379396e-23_rt,  // Ni56
    };
}

namespace Rates
{

    enum NetworkRates
    {
        k_He4_C12_to_O16_reaclib = 1,
        k_He4_O16_to_Ne20_reaclib = 2,
        k_He4_Ne20_to_Mg24_reaclib = 3,
        k_C12_C12_to_He4_Ne20_reaclib = 4,
        k_C12_O16_to_He4_Mg24_reaclib = 5,
        k_O16_O16_to_He4_Si28_reaclib = 6,
        k_C12_Ne20_to_He4_Si28_reaclib = 7,
        k_He4_He4_He4_to_C12_reaclib = 8,
        k_He4_Mg24_to_Si28_removed = 9,
        k_p_Al27_to_Si28_removed = 10,
        k_p_Al27_to_He4_Mg24_removed = 11,
        k_He4_Si28_to_S32_removed = 12,
        k_p_P31_to_S32_removed = 13,
        k_p_P31_to_He4_Si28_removed = 14,
        k_He4_S32_to_Ar36_removed = 15,
        k_p_Cl35_to_Ar36_removed = 16,
        k_p_Cl35_to_He4_S32_removed = 17,
        k_He4_Ar36_to_Ca40_removed = 18,
        k_p_K39_to_Ca40_removed = 19,
        k_p_K39_to_He4_Ar36_removed = 20,
        k_He4_Ca40_to_Ti44_removed = 21,
        k_p_Sc43_to_Ti44_removed = 22,
        k_p_Sc43_to_He4_Ca40_removed = 23,
        k_He4_Ti44_to_Cr48_removed = 24,
        k_He4_Ti44_to_p_V47_removed = 25,
        k_p_V47_to_Cr48_removed = 26,
        k_He4_Cr48_to_Fe52_removed = 27,
        k_He4_Cr48_to_p_Mn51_removed = 28,
        k_p_Mn51_to_Fe52_removed = 29,
        k_He4_Fe52_to_Ni56_removed = 30,
        k_He4_Fe52_to_p_Co55_removed = 31,
        k_p_Co55_to_Ni56_removed = 32,
        k_Mg24_He4_to_Si28_approx = 33,
        k_Si28_to_Mg24_He4_approx = 34,
        k_Si28_He4_to_S32_approx = 35,
        k_S32_to_Si28_He4_approx = 36,
        k_S32_He4_to_Ar36_approx = 37,
        k_Ar36_to_S32_He4_approx = 38,
        k_Ar36_He4_to_Ca40_approx = 39,
        k_Ca40_to_Ar36_He4_approx = 40,
        k_Ca40_He4_to_Ti44_approx = 41,
        k_Ti44_to_Ca40_He4_approx = 42,
        k_Ti44_He4_to_Cr48_approx = 43,
        k_Cr48_to_Ti44_He4_approx = 44,
        k_Cr48_He4_to_Fe52_approx = 45,
        k_Fe52_to_Cr48_He4_approx = 46,
        k_Fe52_He4_to_Ni56_approx = 47,
        k_Ni56_to_Fe52_He4_approx = 48,
        k_O16_to_He4_C12_derived = 49,
        k_Ne20_to_He4_O16_derived = 50,
        k_Mg24_to_He4_Ne20_derived = 51,
        k_C12_to_He4_He4_He4_derived = 52,
        k_He4_Ne20_to_C12_C12_derived = 53,
        k_He4_Mg24_to_C12_O16_derived = 54,
        k_He4_Si28_to_C12_Ne20_derived = 55,
        k_He4_Si28_to_O16_O16_derived = 56,
        k_He4_Mg24_to_p_Al27_removed = 57,
        k_Si28_to_He4_Mg24_removed = 58,
        k_Si28_to_p_Al27_removed = 59,
        k_He4_Si28_to_p_P31_removed = 60,
        k_S32_to_He4_Si28_removed = 61,
        k_S32_to_p_P31_removed = 62,
        k_He4_S32_to_p_Cl35_removed = 63,
        k_Ar36_to_He4_S32_removed = 64,
        k_Ar36_to_p_Cl35_removed = 65,
        k_He4_Ar36_to_p_K39_removed = 66,
        k_Ca40_to_He4_Ar36_removed = 67,
        k_Ca40_to_p_K39_removed = 68,
        k_He4_Ca40_to_p_Sc43_removed = 69,
        k_Ti44_to_He4_Ca40_removed = 70,
        k_Ti44_to_p_Sc43_removed = 71,
        k_Cr48_to_He4_Ti44_removed = 72,
        k_Cr48_to_p_V47_removed = 73,
        k_p_V47_to_He4_Ti44_removed = 74,
        k_Fe52_to_He4_Cr48_removed = 75,
        k_Fe52_to_p_Mn51_removed = 76,
        k_p_Mn51_to_He4_Cr48_removed = 77,
        k_Ni56_to_He4_Fe52_removed = 78,
        k_Ni56_to_p_Co55_removed = 79,
        k_p_Co55_to_He4_Fe52_removed = 80,
        NumRates = k_p_Co55_to_He4_Fe52_removed
    };

    enum ScreenPairs
    {
        k_He4_O16 = 1,
        k_p_K39 = 2,
        k_He4_Ca40 = 3,
        k_He4_Ar36 = 4,
        k_p_Sc43 = 5,
        k_He4_C12 = 6,
        k_He4_Ti44 = 7,
        k_He4_Mg24 = 8,
        k_C12_O16 = 9,
        k_p_V47 = 10,
        k_He4_Cr48 = 11,
        k_C12_C12 = 12,
        k_p_Co55 = 13,
        k_He4_Ne20 = 14,
        k_p_Al27 = 15,
        k_p_Mn51 = 16,
        k_He4_Fe52 = 17,
        k_He4_S32 = 18,
        k_p_P31 = 19,
        k_He4_He4 = 20,
        k_He4_Si28 = 21,
        k_C12_Ne20 = 22,
        k_p_Cl35 = 23,
        k_O16_O16 = 24,
        k_He4_Be8 = 25,
        NumScreenPairs = k_He4_Be8
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "He4_C12_to_O16_reaclib",  // 1,
        "He4_O16_to_Ne20_reaclib",  // 2,
        "He4_Ne20_to_Mg24_reaclib",  // 3,
        "C12_C12_to_He4_Ne20_reaclib",  // 4,
        "C12_O16_to_He4_Mg24_reaclib",  // 5,
        "O16_O16_to_He4_Si28_reaclib",  // 6,
        "C12_Ne20_to_He4_Si28_reaclib",  // 7,
        "He4_He4_He4_to_C12_reaclib",  // 8,
        "He4_Mg24_to_Si28_removed",  // 9,
        "p_Al27_to_Si28_removed",  // 10,
        "p_Al27_to_He4_Mg24_removed",  // 11,
        "He4_Si28_to_S32_removed",  // 12,
        "p_P31_to_S32_removed",  // 13,
        "p_P31_to_He4_Si28_removed",  // 14,
        "He4_S32_to_Ar36_removed",  // 15,
        "p_Cl35_to_Ar36_removed",  // 16,
        "p_Cl35_to_He4_S32_removed",  // 17,
        "He4_Ar36_to_Ca40_removed",  // 18,
        "p_K39_to_Ca40_removed",  // 19,
        "p_K39_to_He4_Ar36_removed",  // 20,
        "He4_Ca40_to_Ti44_removed",  // 21,
        "p_Sc43_to_Ti44_removed",  // 22,
        "p_Sc43_to_He4_Ca40_removed",  // 23,
        "He4_Ti44_to_Cr48_removed",  // 24,
        "He4_Ti44_to_p_V47_removed",  // 25,
        "p_V47_to_Cr48_removed",  // 26,
        "He4_Cr48_to_Fe52_removed",  // 27,
        "He4_Cr48_to_p_Mn51_removed",  // 28,
        "p_Mn51_to_Fe52_removed",  // 29,
        "He4_Fe52_to_Ni56_removed",  // 30,
        "He4_Fe52_to_p_Co55_removed",  // 31,
        "p_Co55_to_Ni56_removed",  // 32,
        "Mg24_He4_to_Si28_approx",  // 33,
        "Si28_to_Mg24_He4_approx",  // 34,
        "Si28_He4_to_S32_approx",  // 35,
        "S32_to_Si28_He4_approx",  // 36,
        "S32_He4_to_Ar36_approx",  // 37,
        "Ar36_to_S32_He4_approx",  // 38,
        "Ar36_He4_to_Ca40_approx",  // 39,
        "Ca40_to_Ar36_He4_approx",  // 40,
        "Ca40_He4_to_Ti44_approx",  // 41,
        "Ti44_to_Ca40_He4_approx",  // 42,
        "Ti44_He4_to_Cr48_approx",  // 43,
        "Cr48_to_Ti44_He4_approx",  // 44,
        "Cr48_He4_to_Fe52_approx",  // 45,
        "Fe52_to_Cr48_He4_approx",  // 46,
        "Fe52_He4_to_Ni56_approx",  // 47,
        "Ni56_to_Fe52_He4_approx",  // 48,
        "O16_to_He4_C12_derived",  // 49,
        "Ne20_to_He4_O16_derived",  // 50,
        "Mg24_to_He4_Ne20_derived",  // 51,
        "C12_to_He4_He4_He4_derived",  // 52,
        "He4_Ne20_to_C12_C12_derived",  // 53,
        "He4_Mg24_to_C12_O16_derived",  // 54,
        "He4_Si28_to_C12_Ne20_derived",  // 55,
        "He4_Si28_to_O16_O16_derived",  // 56,
        "He4_Mg24_to_p_Al27_removed",  // 57,
        "Si28_to_He4_Mg24_removed",  // 58,
        "Si28_to_p_Al27_removed",  // 59,
        "He4_Si28_to_p_P31_removed",  // 60,
        "S32_to_He4_Si28_removed",  // 61,
        "S32_to_p_P31_removed",  // 62,
        "He4_S32_to_p_Cl35_removed",  // 63,
        "Ar36_to_He4_S32_removed",  // 64,
        "Ar36_to_p_Cl35_removed",  // 65,
        "He4_Ar36_to_p_K39_removed",  // 66,
        "Ca40_to_He4_Ar36_removed",  // 67,
        "Ca40_to_p_K39_removed",  // 68,
        "He4_Ca40_to_p_Sc43_removed",  // 69,
        "Ti44_to_He4_Ca40_removed",  // 70,
        "Ti44_to_p_Sc43_removed",  // 71,
        "Cr48_to_He4_Ti44_removed",  // 72,
        "Cr48_to_p_V47_removed",  // 73,
        "p_V47_to_He4_Ti44_removed",  // 74,
        "Fe52_to_He4_Cr48_removed",  // 75,
        "Fe52_to_p_Mn51_removed",  // 76,
        "p_Mn51_to_He4_Cr48_removed",  // 77,
        "Ni56_to_He4_Fe52_removed",  // 78,
        "Ni56_to_p_Co55_removed",  // 79,
        "p_Co55_to_He4_Fe52_removed"  // 80,
    };

}

#endif
