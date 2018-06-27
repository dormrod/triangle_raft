#include "nodes.h"

//##### NODE #####
Node::Node() {
    //default constructor
    id=-1;
}

Node::Node(int idValue, int maxAtomsM, int maxAtomsX, int maxUnits, int maxRings) {
    //constructor
    id=idValue;
    atomsM=Connector(maxAtomsM);
    atomsX=Connector(maxAtomsX);
    units=Connector(maxUnits);
    rings=Connector(maxRings);
}

int Node::addAtomMCnx(int id) {
    //add connection, return 1 if error 0 otherwise
    int status=atomsM.add(id);
    return status;
}

int Node::addAtomXCnx(int id) {
    //add connection, return 1 if error 0 otherwise
    int status=atomsX.add(id);
    return status;
}

int Node::addUnitCnx(int id) {
    //add connection, return 1 if error 0 otherwise
    int status=units.add(id);
    return status;
}

int Node::addRingCnx(int id) {
    //add connection, return 1 if error 0 otherwise
    int status=rings.add(id);
    return status;
}

