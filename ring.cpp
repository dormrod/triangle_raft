#include "ring.h"

//#####  RING #####
Ring::Ring():Node() {} //default constructor

Ring::Ring(int idValue, int maxAtomMCnxs, int maxAtomXCnxs, int maxUnitCnxs, int maxRingCnxs):
        Node(idValue,maxAtomMCnxs,maxAtomXCnxs,maxUnitCnxs,maxRingCnxs) {} //constructor
