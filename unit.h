#ifndef MX2_UNIT_H
#define MX2_UNIT_H

#include <iostream>
#include "connector.h"

using namespace std;

struct Unit {
    //contains connectivity information to atoms, units and rings

    int id; //unique identifier
    int atomM; //metal atom
    Connector atomsX, units, rings; //connections
    bool flag; //generic flag

    Unit(); //constructor
    Unit(int idValue, int maxX, int maxU, int maxR); //constructor

    void setAtomM(int m);
};


#endif //MX2_UNIT_H
