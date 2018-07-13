#ifndef MX2_NETWORK3DS_H
#define MX2_NETWORK3DS_H

#include <iostream>
#include <string>
#include "logfile.h"
#include "crd.h"
#include "networkBase.h"
#include "geom_opt_algs.h"
#include "potentials.h"

using namespace std;

class NetworkCart3DS:public Network<Cart3D> {
    //network class constrained to sphere

protected:
    //Geometry Optimisation
    SteepestDescentArmijo<HLJC3S> optimiser;

    //Virtual Methods To Define
    vector<double> getCrds() override; //get all atom coordinates
    vector<double> getCrds(map<int,int> &globalAtomMap, int n); //get local atom coordinates
    void setCrds(vector<double> &crds) override; //set all atom coordinates
    void setCrds(map<int,int> &globalAtomMap, vector<double> &crds); //set all atom coordinates

public:
    //Constructors
    NetworkCart3DS();
    NetworkCart3DS(string prefix, Logfile &logfile, double additionalParams=0.0); //load network from files

    //Virtual Methods To Define
    void setGO(int it, double ls, double conv, int loc) override; //set up optimiser
    void geometryOptimiseGlobal(vector<double> &potentialModel) override; //optimise geometry with steepest descent
    void geometryOptimiseLocal(vector<double> &potentialModel) override; //optimise geometry with steepest descent
    void writeNetwork(string prefix, Logfile &logfile) override; //write network to files

    int getActiveUnit(string shape, double size) override; //find active unit within shape
    void buildRing(int ringSize, vector<int> &unitPath, vector<double> &potentialModel) override; //build a ring of given size
    void popRing(int ringSize, vector<int> &unitPath) override; //remove last built ring
    void checkOverlap() override; //check for overlap
    bool checkGrowth() override; //check to continue growth
};


#endif //MX2_NETWORK3DS_H
