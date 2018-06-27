#ifndef MX2_NETWORKBASE_H
#define MX2_NETWORKBASE_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "atom.h"
#include "unit.h"
#include "ring.h"

using namespace std;

template <typename CrdT>
class Network {
    //base class for network containing all atom, unit and ring information

protected:
    //key variables
    int nAtoms, nUnits, nRings; //number of atoms, units and rings
    vector< Atom<CrdT> > atoms; //atoms in network (both m and x)
    vector<Unit> units; //triangles in network
    vector<Ring> rings; //rings in network
    map<int,int> atomMap, unitMap, ringMap; //relate ids to references in array

public:
    //constructors
    Network();

    //methods
    void addAtom(Atom<CrdT> atom);
    void addUnit(Unit unit);
    void addRing(Ring ring);
    int addUnitAtomXCnx(int uId, int aId);
    int addUnitRingCnx(int uId, int rId);
    int addUnitUnitCnx(int uId1, int uId2);
    int addRingRingCnx(int rId1, int rId2);

    //write out
    virtual void write(string prefix, Logfile &logfile)=0; //write network to files
};

#include "networkBase.tpp"

#endif //MX2_NETWORKBASE_H
