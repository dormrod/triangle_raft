#include "atom.h"

template <typename CrdT>
Atom<CrdT>::Atom() {
    //default constructor
    id=-1;
}

template <typename CrdT>
Atom<CrdT>::Atom(int idValue, int elem, int cnd) {
    //constructor
    id=idValue;
    element=elem;
    coordination=cnd;
}

template <typename CrdT>
void Atom<CrdT>::setCoordinate(CrdT crd) {
    //set coordinate with given type
    coordinate=crd;
}



