#ifndef MX2_UNIT_H
#define MX2_UNIT_H

#include <iostream>
#include "nodes.h"

using namespace std;

struct Unit: public Node {

    Unit(); //constructor
    Unit(int idValue, int maxAtomMCnxs, int maxAtomXCnxs, int maxUnitCnxs, int maxRingCnxs); //constructor
};


#endif //MX2_UNIT_H
