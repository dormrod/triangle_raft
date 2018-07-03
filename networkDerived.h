#ifndef MX2_NETWORKDERIVED_H
#define MX2_NETWORKDERIVED_H

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
    SteepestDescent<HC2> optimiser;

    //Virtual Methods To Define
    vector<double> getCrds() override; //get all atom coordinates
    void setCrds(vector<double> &crds) override; //set all atom coordinates

public:
    //Constructors
    NetworkCart2D();
    NetworkCart2D(string prefix, Logfile &logfile); //load network from files

    //Virtual Methods To Define
    void setGO(int it, double ls, double conv) override; //set up optimiser
    void geometryOptimise(vector<double> &potentialModel) override; //optimise geometry with steepest descent
    int getActiveUnit(string shape, double size) override; //find active unit within shape
    void buildRing(int ringSize, vector<int> &unitPath, vector<double> &potentialModel) override; //build a ring of given size
    void write(string prefix, Logfile &logfile) override; //write network to files
};


#endif //MX2_NETWORKDERIVED_H
