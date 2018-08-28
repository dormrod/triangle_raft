#ifndef MX2_NETWORKBASE_H
#define MX2_NETWORKBASE_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include "easyIO.h"
#include "logfile.h"
#include "atom.h"
#include "unit.h"
#include "ring.h"
#include "analysis_tools.h"

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
    //Structural
    vector<int> boundaryUnits;
    vector<int> boundaryStatus;
    //Geometry Optimisation
    int optIterations; //number of optimisation iterations
    double defLineInc; //default line search increment
    double energy; //potential energy
    int localExtent, nLocalAtoms; //size of local region, number of atoms in local region
    map<int,int> localAtomMap, globalAtomMap; //maps local to global atoms
    vector<int> flexLocalUnits, fixedLocalUnits, fixedLocalAtoms; //units that make up local region, and fixed atoms
    //Analysis
    bool unitOverlap; //check for overlap of units
    DiscreteDistribution ringStatistics; //ring size distribution for entire network
    map<int,DiscreteDistribution> indRingStatistics; //ring size distributions around each individual ring
    vector<double> aboavWeaireParameters; //alpha, mu and rsq
    ContinuousDistribution bondLenDistXX, bondLenDistMX; //bond length distributions
    map<int,DiscreteDistribution> clusterDistributions; //cluster size distributions for ring sizes
    map<int,bool> percolation; //logs if spanning cluster present for different ring sizes
    vector< col_vector<int> > ringColours; //for visualisation

    //Methods
    virtual vector<double> getCrds()=0; //get all atom coordinates - virtual as have different number of variables and will be faster
    virtual void setCrds(vector<double> &crds)=0; //set all atom coordinates - virtual as have different number of variables and will be faster
    bool checkActiveUnit(int &uId, int sumCheck=12); //checks if active by summing associated atom coordination
    bool checkEdgeUnit(int &uId, int ringCheck=3); //checks if edge by number of associated rings
    void calculateBoundary(); //work out boundary units
    virtual int getActiveUnit(string shape, double size)=0; //find active unit within shape
    void findLocalRegion(int &rId, int nFlexShells); //find local units around given ring

public:
    //Constructors
    Network();

    //Setters
    virtual void setGO(int it, double ls, double conv, int loc)=0; //virtual as set up optimiser with different potential types

    //Getters
    int getNRings();
    double getEnergy();
    int getIterations();

    //Methods
    //Build Network
    void addAtom(Atom<CrdT> atom);
    void addUnit(Unit unit);
    void addRing(Ring ring);
    void delAtom();
    void delUnit();
    void delRing();
    int addUnitAtomXCnx(int uId, int aId);
    int addUnitRingCnx(int uId, int rId);
    int addUnitUnitCnx(int uId1, int uId2);
    int addRingRingCnx(int rId1, int rId2);
    void delUnitUnitCnx(int uId1, int uId2);
    void delUnitRingCnx(int uId, int rId);
    void delRingRingCnx(int rId1, int rId2);
    void changeUnitAtomXCnx(int uId, int aId1, int aId2);
    bool trialRing(int ringSize, vector<int> &unitPath, vector<double> &potentialModel); //test a trial ring of given size
    void acceptRing(int ringSize, vector<int> &unitPath, vector<double> &potentialModel); //accept a ring of given size
    void clean(); //clean network of dead atoms
    virtual void buildRing(int ringSize, vector<int> &unitPath, vector<double> &potentialModel)=0; //build a ring of given size
    virtual void buildRing0(vector<int> &unitPath)=0; //build a ring of same size as unit path
    virtual void popRing(int ringSize, vector<int> &unitPath)=0; //remove last built ring
    virtual void popRing0(vector<int> &unitPath)=0; //remove last built ring of same size as unit path
    virtual bool checkGrowth()=0; //check to continue growth
    virtual bool checkLocalGrowth(int rId)=0; //check if acceptable local growth

    //Search Network
    vector<int> getBoundarySection(int startId, bool direction); //find section of unit boundary in given direction
    //Optimise Network
    virtual void geometryOptimiseGlobal(vector<double> &potentialModel)=0;
    virtual void geometryOptimiseLocal(vector<double> &potentialModel)=0;
    //Analyse Network
    void calculateRingStatistics(); //ring stats analysis
    void calculateBondDistributions(); //bond len/angle distributions
    virtual void checkOverlap()=0; //check for overlap
    virtual void calculatePercolation(string shape)=0; //clusters and percolation

    //Write Out
    void write(string prefix, bool special, Logfile &logfile); //write out to files
    virtual void writeNetwork(string prefix, Logfile &logfile)=0; //write network to files
    virtual void writeNetworkSpecial(string prefix, Logfile &logfile)=0; //write network in different format
    void writeAnalysis(string prefix, Logfile &logfile); //write analysis out to files
    void kill(string prefix, Logfile &logfile); //prepare network for early termination and write out
};

#include "networkBase.tpp"

#endif //MX2_NETWORKBASE_H
