#include "netBaseUnits.h"

//##### ATOM #####
template <typename CrdT>
Atom<CrdT>::Atom() {
    //default constructor
    id=-1;
}

template <typename CrdT>
Atom<CrdT>::Atom(int initId, int initElement) {
    //constructor
    id=initId;
    element=initElement;
}

template <typename CrdT>
Atom<CrdT>::~Atom() {
    //destructor
}

template <typename CrdT>
Atom<CrdT>::Atom(const Atom &source) {
    //copy constructor

    //shallow copies
    id=source.id;
    element=source.element;
}

template <typename CrdT>
Atom<CrdT>& Atom<CrdT>::operator=(const Atom<CrdT> &source) {
    //overload assignment operator
    if (this == &source) return *this;

    //shallow copies
    id=source.id;
    element=source.element;

    return *this;
}

//##### CRYSTAL UNIT #####

template <typename CrdT>
CrystalUnit<CrdT>::CrystalUnit() {
    //default constructor
    id=-1;
}

template <typename CrdT>
CrystalUnit<CrdT>::CrystalUnit(int initId, int initNX) {
    //constructor
    id=initId;
    nX=initNX;
    mAtom=-1;
    xAtoms=new int[nX]();
}

template <typename CrdT>
CrystalUnit<CrdT>::~CrystalUnit() {
    //destructor
    delete[] xAtoms;
}

template <typename CrdT>
CrystalUnit<CrdT>::CrystalUnit(const CrystalUnit &source) {
    //copy constructor

    //shallow copies
    id=source.id;
    nX=source.nX;
    mAtom=source.mAtom;

    //deep copies
    xAtoms=new int[nX]();
    for(int i=0; i<nX; ++i) xAtoms[i]=source.xAtoms[i];
}

template <typename CrdT>
CrystalUnit<CrdT>& CrystalUnit<CrdT>::operator=(const CrystalUnit<CrdT> &source) {
    //overload assigment operator
    if (this == &source) return *this;

    //shallow copies
    id=source.id;
    nX=source.nX;
    mAtom=source.mAtom;

    //deep copies
    xAtoms=new int[nX]();
    for(int i=0; i<nX; ++i) xAtoms[i]=source.xAtoms[i];

    return *this;
}

//##### RING #####
template <typename CrdT>
Ring<CrdT>::Ring() {
    //default constructor
    id=-1;
}

template <typename CrdT>
Ring<CrdT>::Ring(int initId, int initN) {
    //constructor
    id=initId;
    n=initN;
    units=new int[n]();
}

template <typename CrdT>
Ring<CrdT>::~Ring() {
    //destructor
    delete[] units;
}

template <typename CrdT>
Ring<CrdT>::Ring(const Ring &source) {
    //copy constructor

    //shallow copies
    id=source.id;
    n=source.n;

    //deep copies
    units=new int[n]();
    for(int i=0; i<n; ++i) units[i]=source.units[i];
}

template <typename CrdT>
Ring<CrdT>& Ring<CrdT>::operator=(const Ring<CrdT> &source) {
    //overload assignment operator
    if (this == &source) return *this;

    //shallow copies
    id=source.id;
    n=source.n;

    //deep copies
    units=new int[n]();
    for(int i=0; i<n; ++i) units[i]=source.units[i];

    return *this;
}
