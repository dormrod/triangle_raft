#ifndef MX2_ATOM_H
#define MX2_ATOM_H

#include <iostream>
#include "nodes.h"

template <typename CrdT>
struct Atom: public Node {
    //contains element and coordinate information of given type

    int element; //atomic number
    CrdT coordinate; //of given type

    Atom(); //constructor
    Atom(int elem, int idValue, int maxAtomMCnxs, int maxAtomXCnxs, int maxUnitCnxs, int maxRingCnxs); //construct with element
    Atom(int elem, CrdT crd, int idValue, int maxAtomMCnxs, int maxAtomXCnxs, int maxUnitCnxs, int maxRingCnxs); //construct with element and coordinate

    void setCoordinate(CrdT crd); //set coordinate
};

#include "atom.tpp"

#endif //MX2_ATOM_H
