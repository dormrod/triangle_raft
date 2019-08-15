#ifndef MX2_ATOM_H
#define MX2_ATOM_H

#include <iostream>

template <typename CrdT>
struct Atom {
    //contains element, coordinate and unit connectivity information

    int id; //identifier
    int element; //atomic no.
    int coordination; //number of bonds
    CrdT coordinate; //of give type

    Atom(); //constructor
    Atom(int idValue, int elem, int cnd); //construct with id, element and coordination

    void setCoordinate(CrdT crd); //set coordinate
};

#include "atom.tpp"

#endif //MX2_ATOM_H
