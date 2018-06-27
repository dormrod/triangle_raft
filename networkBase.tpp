#include "networkBase.h"


//##### NETWORK BASE #####

template <typename CrdT>
Network<CrdT>::Network() {
    //default constructor
    nAtoms=0;
    nGeomUnits=0;
    nRings=0;
    atoms.clear();
    geomUnits.clear();
    rings.clear();
}

