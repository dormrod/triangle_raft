#ifndef MX2_RING_H
#define MX2_RING_H

#include <iostream>
#include "connector.h"

using namespace std;

struct Ring{
    //contains unit and ring connection information

    int id;//unique identifier
    Connector units, rings; //connections

    Ring(); //constructor
    Ring(int idValue, int maxU, int maxR); //constructor
};


#endif //MX2_RING_H
