#ifndef NETWORK_PROPERTIES_H
#define NETWORK_PROPERTIES_H
#include <vector>
#include <string>
#include <amrex_bridge.H>

constexpr int NumSpec = 13;

constexpr Real aion[NumSpec] = {
       4.0, // 0
      12.0, // 1
      16.0, // 2
      20.0, // 3
      24.0, // 4
      28.0, // 5
      32.0, // 6
      36.0, // 7
      40.0, // 8
      44.0, // 9
      48.0, // 10
      52.0, // 11
      56.0, // 12
 };

constexpr Real aion_inv[NumSpec] = {
    1.0/   4.0, // 0
    1.0/  12.0, // 1
    1.0/  16.0, // 2
    1.0/  20.0, // 3
    1.0/  24.0, // 4
    1.0/  28.0, // 5
    1.0/  32.0, // 6
    1.0/  36.0, // 7
    1.0/  40.0, // 8
    1.0/  44.0, // 9
    1.0/  48.0, // 10
    1.0/  52.0, // 11
    1.0/  56.0, // 12
 };

constexpr Real zion[NumSpec] = {
       2.0, // 0
       6.0, // 1
       8.0, // 2
      10.0, // 3
      12.0, // 4
      14.0, // 5
      16.0, // 6
      18.0, // 7
      20.0, // 8
      22.0, // 9
      24.0, // 10
      26.0, // 11
      28.0, // 12
 };

static const std::vector<std::string> spec_names = {
    "He4", // 0
    "C12", // 1
    "O16", // 2
    "Ne20", // 3
    "Mg24", // 4
    "Si28", // 5
    "S32", // 6
    "Ar36", // 7
    "Ca40", // 8
    "Ti44", // 9
    "Cr48", // 10
    "Fe52", // 11
    "Ni56", // 12
 };

namespace Species {
  enum NetworkSpecies {
    He4=1,
    C12,
    O16,
    Ne20,
    Mg24,
    Si28,
    S32,
    Ar36,
    Ca40,
    Ti44,
    Cr48,
    Fe52,
    Ni56,
  };
}

#endif
