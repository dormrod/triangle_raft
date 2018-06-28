#include "networkBase.h"


//##### NETWORK BASE #####

template <typename CrdT>
Network<CrdT>::Network() {
    //default constructor
    nAtoms=0;
    nUnits=0;
    nRings=0;
    energy=numeric_limits<double>::infinity();
    optIterations=-1;
    atoms.clear();
    units.clear();
    rings.clear();
}

template <typename CrdT>
int Network<CrdT>::getNRings() {
    //return number of rings in network
    return nRings;
}

template <typename CrdT>
void Network<CrdT>::addAtom(Atom<CrdT> atom) {
    //add atom to network and update map
    atoms.push_back(atom);
    ++nAtoms;
}

template <typename CrdT>
void Network<CrdT>::addUnit(Unit unit) {
    //add unit to network and update map
    units.push_back(unit);
    ++nUnits;
}

template <typename CrdT>
void Network<CrdT>::addRing(Ring ring) {
    //add ring to network and update map
    rings.push_back(ring);
    ++nRings;
}

template <typename CrdT>
int Network<CrdT>::addUnitAtomXCnx(int uId, int aId) {
    //add atom x id to unit
    int status=units[uId].atomsX.add(aId);
    return status;
}

template <typename CrdT>
int Network<CrdT>::addUnitRingCnx(int uId, int rId) {
    //add mutual unit-ring connection
    int uStatus=units[uId].rings.add(rId);
    int rStatus=rings[rId].units.add(uId);
    return uStatus+rStatus;
}

template <typename CrdT>
int Network<CrdT>::addUnitUnitCnx(int uId1, int uId2) {
    //add mutual unit-unit connection
    int status1=units[uId1].units.add(uId2);
    int status2=units[uId2].units.add(uId1);
    return status1+status2;
}

template <typename CrdT>
int Network<CrdT>::addRingRingCnx(int rId1, int rId2) {
    //add mutual unit-unit connection
    int status1=rings[rId1].rings.add(rId2);
    int status2=rings[rId2].rings.add(rId1);
    return status1+status2;
}

