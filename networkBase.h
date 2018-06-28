#ifndef MX2_NETWORKBASE_H
#define MX2_NETWORKBASE_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "easyIO.h"
#include "atom.h"
#include "unit.h"
#include "ring.h"

using namespace std;

template <typename CrdT>
class Network {
    //base class for network containing all atom, unit and ring information

protected:
    //Key Variables
    int nAtoms, nUnits, nRings; //number of atoms, units and rings
    vector< Atom<CrdT> > atoms; //atoms in network (both m and x)
    vector<Unit> units; //triangles in network
    vector<Ring> rings; //rings in network

    //Additional Variables
    //Geometry Optimisation
    int optIterations; //number of optimisation iterations
    double energy; //potential energy


    //Methods
    virtual vector<double> getCrds()=0; //get all atom coordinates - virtual as have different number of variables and will be faster
    virtual void setCrds(vector<double> &crds)=0; //set all atom coordinates - virtual as have different number of variables and will be faster

public:
    //Constructors
    Network();

    //Setters
    virtual void setGO(int it, double ls, double conv)=0; //virtual as set up optimiser with different potential types

    //Getters
    int getNRings();

    //Methods
    //Build Network
    void addAtom(Atom<CrdT> atom);
    void addUnit(Unit unit);
    void addRing(Ring ring);
    int addUnitAtomXCnx(int uId, int aId);
    int addUnitRingCnx(int uId, int rId);
    int addUnitUnitCnx(int uId1, int uId2);
    int addRingRingCnx(int rId1, int rId2);
    //Optimise Network
    virtual void geometryOptimise(vector<double> &potentialModel)=0;

    //Write Out
    virtual void write(string prefix, Logfile &logfile)=0; //write network to files
};

#include "networkBase.tpp"

#endif //MX2_NETWORKBASE_H
