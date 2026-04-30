#include <actual_network.H>
#include <amrex_bridge.H>

namespace network
{
    Array1D<Real, 1, NumSpec> bion;
    Array1D<Real, 1, NumSpec> mion;
}

void actual_network_init()
{
    using namespace Species;
    using namespace network;

    // binding energies per nucleon in MeV
    Array1D<Real, 1, NumSpec> ebind_per_nucleon;

    ebind_per_nucleon(H1) = 0.0_rt;
    ebind_per_nucleon(He4) = 7.073915614499924_rt;
    ebind_per_nucleon(C12) = 7.680144581999836_rt;
    ebind_per_nucleon(N13) = 7.238863542153782_rt;
    ebind_per_nucleon(N14) = 7.475614810571252_rt;
    ebind_per_nucleon(O16) = 7.97620721324995_rt;
    ebind_per_nucleon(F18) = 7.6316390264443426_rt;
    ebind_per_nucleon(Ne20) = 8.032241192000038_rt;
    ebind_per_nucleon(Ne21) = 7.971713797142694_rt;
    ebind_per_nucleon(Na22) = 7.915662309272772_rt;
    ebind_per_nucleon(Na23) = 8.111493582782714_rt;
    ebind_per_nucleon(Mg24) = 8.26071033199984_rt;
    ebind_per_nucleon(Al27) = 8.331553230814784_rt;
    ebind_per_nucleon(Si28) = 8.447744478428442_rt;
    ebind_per_nucleon(P31) = 8.481167740645281_rt;
    ebind_per_nucleon(S32) = 8.493130116374914_rt;
    ebind_per_nucleon(Ar36) = 8.51990963755553_rt;
    ebind_per_nucleon(Ca40) = 8.55130463200003_rt;
    ebind_per_nucleon(Ti44) = 8.53352185472729_rt;
    ebind_per_nucleon(Cr48) = 8.572248748666576_rt;
    ebind_per_nucleon(Fe52) = 8.609608043538314_rt;
    ebind_per_nucleon(Ni56) = 8.642780296285569_rt;

    // convert to binding energies per nucleus in MeV
    for (int i = 1; i <= NumSpec; ++i) {
        bion(i) = ebind_per_nucleon(i) * aion[i-1];
    }

    // Set the mass -- this will be in grams
    mion(H1) = 1.6735328377636005e-24_rt;
    mion(He4) = 6.646479071584587e-24_rt;
    mion(C12) = 1.99264687992e-23_rt;
    mion(N13) = 2.1596537049448796e-23_rt;
    mion(N14) = 2.3252651436495096e-23_rt;
    mion(O16) = 2.6560180592333686e-23_rt;
    mion(F18) = 2.989125964092377e-23_rt;
    mion(Ne20) = 3.3198227947612416e-23_rt;
    mion(Ne21) = 3.4861102572650884e-23_rt;
    mion(Na22) = 3.652262279854593e-23_rt;
    mion(Na23) = 3.817541002484691e-23_rt;
    mion(Mg24) = 3.9828098739467446e-23_rt;
    mion(Al27) = 4.480389861070653e-23_rt;
    mion(Si28) = 4.6456779473820677e-23_rt;
    mion(P31) = 5.14331418367544e-23_rt;
    mion(S32) = 5.309087322384128e-23_rt;
    mion(Ar36) = 5.972551377884467e-23_rt;
    mion(Ca40) = 6.635944331004904e-23_rt;
    mion(Ti44) = 7.299678247096977e-23_rt;
    mion(Cr48) = 7.962953983065421e-23_rt;
    mion(Fe52) = 8.626187166893794e-23_rt;
    mion(Ni56) = 9.289408870379396e-23_rt;

}
