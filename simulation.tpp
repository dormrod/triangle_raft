#include "simulation.h"

template <typename CrdT, typename NetT>
Simulation<CrdT,NetT>::Simulation() {
    //default constructor
}

template <typename CrdT, typename NetT>
Simulation<CrdT,NetT>::Simulation(Logfile &logfile) {
    //constructor
    logfile.log("Initialising network parameters","","",0,false);
}

template <typename CrdT, typename NetT>
void Simulation<CrdT,NetT>::setIO(string in, string out, Logfile &logfile) {
    //set io parameters
    prefixIn=in;
    prefixOut=out;
    logfile.log("Initialised: ","IO","",1,false);
}

template <typename CrdT, typename NetT>
void Simulation<CrdT,NetT>::setNP(int targRings, int basicMin, int basicMax, string geom, Logfile &logfile) {
    //set network properties
    nTargetRings=targRings;
    basicMinSize=basicMin;
    basicMaxSize=basicMax;
    if(geom.substr(0,2)=="2D") dimensionality=2;
    else if(geom.substr(0,2)=="3D") dimensionality=3;
    else logfile.errorlog("Geometry code incorrect","critical");
    growthGeometry=geom.substr(2,2);
    logfile.log("Initialised: ","network properties","",1,false);
}

template <typename CrdT, typename NetT>
void Simulation<CrdT,NetT>::setMC(int seed, double temperature, Logfile &logfile) {
    //set up monte carlo generator
    monteCarlo=MonteCarlo(seed,temperature,logfile);
}

template <typename CrdT, typename NetT>
void Simulation<CrdT,NetT>::setPM(double kMX, double r0MX, double kXX, double a0XX, double kMM, double a0MM, Logfile &logfile) {
    //set potential model as M-X, X-X, M-M k and r0 values for harmonics
    potentialModel.resize(6);
    potentialModel[0]=kMX;
    potentialModel[1]=r0MX;
    potentialModel[2]=kXX;
    potentialModel[4]=kMM;
    //convert angles to equivalent distances
    double r0=2.0*kMX*sin(a0XX*0.5*M_PI/180.0);
    potentialModel[3]=r0;
    r0=2.0*kMX*sin(a0MM*0.5*M_PI/180.0);
    potentialModel[5]=r0;
    logfile.log("Initialised: ","potential model","",1,false);
}

template <typename CrdT, typename NetT>
void Simulation<CrdT,NetT>::setGO(bool global0, bool global1, int it, double ls, double conv, Logfile &logfile) {
    //set geometery optimisation parameters
    globalPreGO=global0;
    globalPostGO=global1;
    goMaxIterations=it;
    goLineSeachInc=ls;
    goConvergence=conv;
    logfile.log("Initialised: ","geometry optimisation","",1,false);
}

//template <typename CrdT, typename NetT>
//Simulation<CrdT,NetT>::
//template <typename CrdT, typename NetT>
//Simulation<CrdT,NetT>::
//template <typename CrdT, typename NetT>
//Simulation<CrdT,NetT>::
//template <typename CrdT, typename NetT>
//Simulation<CrdT,NetT>::
//template <typename CrdT, typename NetT>
//Simulation<CrdT,NetT>::
//template <typename CrdT, typename NetT>
//Simulation<CrdT,NetT>::
