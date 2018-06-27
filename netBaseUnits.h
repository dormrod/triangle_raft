#ifndef MX2_CRD_H
#define MX2_CRD_H

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

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

    //overloaded operators
    Atom<CrdT>& operator=(const Atom<CrdT> &source);
};

template <typename CrdT>
struct CrystalUnit{
    //base unit of crystal e.g. triangle in 2D, tetrahedron in 3D

    //key variables
    int id; //specific identifier
    int nX; //number of x atoms
    int mAtom; //id of metal atom
    int *xAtoms; //ids of x atoms

    //constructors, destructors
    CrystalUnit();
    CrystalUnit(int initId, int initNX);
    CrystalUnit(const CrystalUnit &source);
    ~CrystalUnit();

    //overloaded operators
    CrystalUnit<CrdT>& operator=(const CrystalUnit<CrdT> &source);
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
    //overloaded operators
    Ring<CrdT>& operator=(const Ring<CrdT> &source);
};

#include "netBaseUnits.tpp"

#endif //MX2_CRD_H
