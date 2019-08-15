#ifndef MX2_NETWORK2DC_H
#define MX2_NETWORK2DC_H

#include <iostream>
#include <string>
#include "logfile.h"
#include "crd.h"
#include "networkBase.h"
#include "geom_opt_algs.h"
#include "potentials.h"

using namespace std;

class NetworkCart2D: public Network<Cart2D> {
    //network class using two dimensional cartesian coordinate
protected:
    //Geometry Optimisation
    SteepestDescentArmijo<HLJC2> optimiser;

    //Virtual Methods To Define
    vector<double> getCrds() override; //get all atom coordinates
    vector<double> getCrds(map<int,int> &globalAtomMap, int n); //get local atom coordinates
    void setCrds(vector<double> &crds) override; //set all atom coordinates
    void setCrds(map<int,int> &globalAtomMap, vector<double> &crds); //set all atom coordinates

public:
    //Constructors
    NetworkCart2D();
    NetworkCart2D(string prefix, Logfile &logfile, double additionalParams=0.0); //load network from files

    //Virtual Methods To Define
    void setGO(int it, double ls, double conv, int loc) override; //set up optimiser
    void geometryOptimiseGlobal(vector<double> &potentialModel) override; //optimise geometry with steepest descent
    void geometryOptimiseLocal(vector<double> &potentialModel) override; //optimise geometry with steepest descent
    int getActiveUnit(string shape, double size) override; //find active unit within shape
    void buildRing(int ringSize, vector<int> &unitPath, vector<double> &potentialModel) override; //build a ring of given size
    void buildRing0(vector<int> &unitPath) override; //build a ring of same size as unit path
    void popRing(int ringSize, vector<int> &unitPath) override; //remove last built ring
    void popRing0(vector<int> &unitPath) override; //remove last built ring of same size as unit path
    void writeNetwork(string prefix, Logfile &logfile) override; //write network to files
    void writeNetworkSpecial(string prefix, Logfile &logfile) override; //write network in different format
    void checkOverlap() override; //check for overlap
    bool checkGrowth() override; //check to continue growth
    bool checkLocalGrowth(int rId) override; //check if acceptable local growth
    void calculatePercolation(string shape) override; //clusters and percolation
};

#endif //MX2_NETWORK2DC_H
