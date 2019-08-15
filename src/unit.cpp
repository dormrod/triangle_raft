#include "unit.h"

//#####  UNIT #####
Unit::Unit(){
    //default constructor
    id=-1;
}

Unit::Unit(int idValue, int maxX, int maxU, int maxR) {
    //constructor
    id=idValue;
    atomsX=Connector(maxX);
    units=Connector(maxU);
    rings=Connector(maxR);
    flag=false;
}

void Unit::setAtomM(int m) {
    //set id of metal atom
    atomM=m;
}

