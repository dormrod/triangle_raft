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
    coordinate=source.coordinate;

}

template <typename CrdT>
Atom<CrdT>& Atom<CrdT>::operator=(const Atom<CrdT> &source) {
    //overload assignment operator
    if (this == &source) return *this;

    //shallow copies
    id=source.id;
    element=source.element;
    coordinate=source.coordinate;

    return *this;
}

template <typename CrdT>
void Atom<CrdT>::setCoordinate(CrdT c) {
    //set coordinate using given type
    coordinate=c;
}

//##### GEOMETRICAL UNIT #####

template <typename CrdT>
GeometricalUnit<CrdT>::GeometricalUnit() {
    //default constructor
    id=-1;
}

template <typename CrdT>
GeometricalUnit<CrdT>::GeometricalUnit(int initId, int initNX) {
    //constructor
    id=initId;
    nX=initNX;
    mAtom=-1;
    xAtoms=new int[nX]();
}

template <typename CrdT>
GeometricalUnit<CrdT>::~GeometricalUnit() {
    //destructor
    delete[] xAtoms;
}

template <typename CrdT>
GeometricalUnit<CrdT>::GeometricalUnit(const GeometricalUnit &source) {
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
GeometricalUnit<CrdT>& GeometricalUnit<CrdT>::operator=(const GeometricalUnit<CrdT> &source) {
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

template <typename CrdT>
void GeometricalUnit<CrdT>::setAtomM(int m) {
    //set m atom id
    mAtom=m;
}

template <typename CrdT>
void GeometricalUnit<CrdT>::setAtomsX(vector<int> x) {
    //set x atoms from vector
    for(int i=0; i<nX; ++i) xAtoms[i]=x[i];
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

template <typename CrdT>
void Ring<CrdT>::setUnits(vector<int> u) {
    //set units from vector
    for(int i=0; i<n; ++i) units[i]=u[i];
}
