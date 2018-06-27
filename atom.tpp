#include "atom.h"

template <typename CrdT>
Atom<CrdT>::Atom():Node() {} //default constructor

template <typename CrdT>
Atom<CrdT>::Atom(int elem, int idValue, int maxAtomMCnxs, int maxAtomXCnxs, int maxUnitCnxs, int maxRingCnxs):
        Node(idValue,maxAtomMCnxs,maxAtomXCnxs,maxUnitCnxs,maxRingCnxs){
    //constructor
    element=elem;
}

template <typename CrdT>
Atom<CrdT>::Atom(int elem, CrdT crd, int idValue, int maxAtomMCnxs, int maxAtomXCnxs, int maxUnitCnxs, int maxRingCnxs):
        Node(idValue,maxAtomMCnxs,maxAtomXCnxs,maxUnitCnxs,maxRingCnxs) {
    //constructor
    element=elem;
    coordinate=crd;
}

template <typename CrdT>
void Atom<CrdT>::setCoordinate(CrdT crd) {
    //set coordinate with given type
    coordinate=crd;
}



