#ifndef MX2_NODES_H
#define MX2_NODES_H

#include <iostream>
#include "connector.h"

using namespace std;

struct Node {
    //contains connections to atoms, units and rings

    //key variables
    int id; //unique identifier
    Connector atomsM, atomsX, units, rings; //connections

    //constructor
    Node();
    Node(int idValue, int maxAtomsM, int maxAtomsX, int maxUnits, int maxRings);

    //methods
    int addAtomMCnx(int id);
    int addAtomXCnx(int id);
    int addUnitCnx(int id);
    int addRingCnx(int id);
};


#endif //MX2_NODES_H
