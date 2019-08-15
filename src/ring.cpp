#include "ring.h"

//#####  RING #####
Ring::Ring(){
    //default constructor
    id=-1;
}

Ring::Ring(int idValue, int maxU, int maxR) {
    //constructor
    id=idValue;
    units=Connector(maxU);
    rings=Connector(maxR);
}
