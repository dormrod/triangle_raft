#include "unit.h"

//#####  UNIT #####
Unit::Unit():Node() {} //default constructor

Unit::Unit(int idValue, int maxAtomMCnxs, int maxAtomXCnxs, int maxUnitCnxs, int maxRingCnxs):
        Node(idValue,maxAtomMCnxs,maxAtomXCnxs,maxUnitCnxs,maxRingCnxs) {} //constructor