#ifndef MX2_NETWORKBASE_H
#define MX2_NETWORKBASE_H

#include <iostream>
#include <vector>
#include <string>
#include "netBaseUnits.h"

using namespace std;

template <typename CrdT>
class Network {
    //base class for network containing all atom, crystal unit and ring information

protected:
    //key variables
    int nAtoms, nGeomUnits, nRings; //number of atoms, geom units and rings
    vector< Atom<CrdT> > atoms; //atoms in network (both m and x)
    vector< GeometricalUnit<CrdT> > geomUnits; //triangles in network
    vector< Ring<CrdT> > rings; //rings in network

public:
    //constructors
    Network();

    //write out
    virtual void write(string prefix, Logfile &logfile)=0; //write network to files
};

#include "networkBase.tpp"

#endif //MX2_NETWORKBASE_H
