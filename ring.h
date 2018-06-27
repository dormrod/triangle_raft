#ifndef MX2_RING_H
#define MX2_RING_H

#include <iostream>
#include "nodes.h"

using namespace std;

struct Ring: public Node {

    Ring(); //constructor
    Ring(int idValue, int maxAtomMCnxs, int maxAtomXCnxs, int maxUnitCnxs, int maxRingCnxs); //constructor
};


#endif //MX2_RING_H
