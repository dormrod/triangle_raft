#ifndef MX2_NETBASEUNITS_H
#define MX2_NETBASEUNITS_H

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "connector.h"

using namespace std;

template <typename CrdT>
struct Atom{
    //contains positional and element information

    //key variables
    int id; //specific identifier
    int element; //atomic number
    CrdT coordinate; //coordinate of template type

    //constructors, destructors
    Atom();
    Atom(int initId, int initElement);
    Atom(const Atom &source);
    ~Atom();

    //setters
    void setCoordinate(CrdT c);

    //overloaded operators
    Atom<CrdT>& operator=(const Atom<CrdT> &source);
};

template <typename CrdT>
struct GeometricalUnit{
    //base geometry e.g. triangle in 2D, tetrahedron in 3D

    //key variables
    int id; //specific identifier
    int nX; //number of x atoms
    int mAtom; //id of metal atom
    int *xAtoms; //ids of x atoms

    //constructors, destructors
    GeometricalUnit();
    GeometricalUnit(int initId, int initNX);
    GeometricalUnit(const GeometricalUnit &source);
    ~GeometricalUnit();

    //setters
    void setAtomM(int m);
    void setAtomsX(vector<int> x);

    //overloaded operators
    GeometricalUnit<CrdT>& operator=(const GeometricalUnit<CrdT> &source);
};

template <typename CrdT>
struct Ring{
    //rings formed from crystal units

    //key variables
    int id; //specific identifier
    int n; //size of ring
    int *units; //ids of crystal units that form ring

    //constructors, destructors
    Ring();
    Ring(int initId, int initN);
    Ring(const Ring &source);
    ~Ring();

    //setters
    void setUnits(vector<int> u);

    //overloaded operators
    Ring<CrdT>& operator=(const Ring<CrdT> &source);
};

#include "netBaseUnits.tpp"

#endif //MX2_NETBASEUNITS_H
