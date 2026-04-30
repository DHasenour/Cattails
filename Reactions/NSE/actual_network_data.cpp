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

    ebind_per_nucleon(He4) = 7.073915614499924_rt;
    ebind_per_nucleon(C12) = 7.680144581999836_rt;
    ebind_per_nucleon(O16) = 7.97620721324995_rt;
    ebind_per_nucleon(Ne20) = 8.032241192000038_rt;
    ebind_per_nucleon(Mg24) = 8.26071033199984_rt;
    ebind_per_nucleon(Si28) = 8.447744478428442_rt;
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
    for (int i = 1; i <= NumSpec; ++i) {
        mion(i) = (aion[i-1] - zion[i-1]) * C::m_n + zion[i-1] * (C::m_p + C::m_e) - bion(i) * C::MeV2gr;
    }

}
