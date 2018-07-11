#include "simulation.h"

//##### INITIALISATION #####
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
    nBasicRingSizes=basicMax-basicMin+1;
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
void Simulation<CrdT,NetT>::setPM(double kMX, double r0MX, double kXX, double a0XX, double kMM, double a0MM, double kLJ, double r0LJ, Logfile &logfile) {
    //set potential model as M-X, X-X, M-M k and r0 values for harmonics
    potentialModel.resize(8);
    potentialModel[0]=kMX;
    potentialModel[1]=r0MX;
    potentialModel[2]=kXX;
    potentialModel[4]=kMM;
    //convert angles to equivalent distances
    double r0=2.0*kMX*sin(a0XX*0.5*M_PI/180.0);
    potentialModel[3]=r0;
    r0=2.0*kMX*sin(a0MM*0.5*M_PI/180.0);
    potentialModel[5]=r0;
    potentialModel[6]=kLJ;
    potentialModel[7]=r0LJ;
//    potentialModel[8]=2.0*sqrt(pow(r0MX,2)-pow(potentialModel[3]*0.5,2));
    logfile.log("Initialised: ","potential model","",1,false);
}

template <typename CrdT, typename NetT>
void Simulation<CrdT,NetT>::setGO(bool global0, bool global1, int it, double ls, double conv, int loc, Logfile &logfile) {
    //set geometery optimisation parameters
    globalPreGO=global0;
    globalPostGO=global1;
    goMaxIterations=it;
    goLineSeachInc=ls;
    goConvergence=conv;
    goLocalExtent=loc;
    logfile.log("Initialised: ","geometry optimisation","",1,false);
}

//##### MAIN #####
template <typename CrdT, typename NetT>
void Simulation<CrdT,NetT>::run(Logfile &logfile) {
    //main control of simulation

    loadNetwork(logfile);

    growNetwork(logfile);

    analyseNetwork(logfile);

    writeNetwork(logfile);
}

//##### LOAD #####
template <typename CrdT, typename NetT>
void Simulation<CrdT,NetT>::loadNetwork(Logfile &logfile) {
    //read in network from file

    logfile.log("Intialisation complete","","",0,true);
    masterNetwork=NetT(prefixIn,logfile);
    masterNetwork.setGO(goMaxIterations,goLineSeachInc,goConvergence,goLocalExtent);
    if(globalPreGO) masterNetwork.geometryOptimiseGlobal(potentialModel);
}

//##### GROW #####
template <typename CrdT, typename NetT>
void Simulation<CrdT,NetT>::growNetwork(Logfile &logfile) {
    //grow network using monte carlo process
    logfile.log("Network growth begun after","","sec",0,false);

    killGrowth=false;
    energyCutoff=100.0*potentialModel[0];
    int nRings=masterNetwork.getNRings();
    if(nRings<nTargetRings){//only if network needs growing
        do{
            int activeUnit = selectActiveUnit();
            vector<int> unitPath = selectUnitPath(activeUnit);
            addBasicRing(unitPath);
            ++nRings;
            if(nRings%100==0){
                cout<<nRings<<endl;
                logfile.log(to_string(nRings)+" rings, time elapsed: ","","sec",1,false);
            }
            cout<<nRings<<endl;
            if(killGrowth){
                masterNetwork.writeNetwork(prefixOut,logfile);
                logfile.errorlog("Growth prematurely killed due to excessive energy","critical");
                break;
            }
        }while(nRings<nTargetRings);
    }
    if(globalPostGO){
        logfile.log("Performing global geometry optimisation","","",1,false);
        masterNetwork.geometryOptimiseGlobal(potentialModel);
    }
    logfile.log("Network growth complete","","",0,true);
}

template <typename CrdT, typename NetT>
int Simulation<CrdT,NetT>::selectActiveUnit() {
    //find unit within geometrical shape which has dangling bonds
    int activeUnit;
    double regionSize=1.0;
    for(;;){
        activeUnit=masterNetwork.getActiveUnit(growthGeometry,regionSize);
        if(activeUnit==-1) regionSize+=1.0;
        else break;
    }
    return activeUnit;
}

template <typename CrdT, typename NetT>
vector<int> Simulation<CrdT,NetT>::selectUnitPath(int activeUnit) {
    //find both paths which could form new ring, choose one due to some criteria
    vector<int> unitPath, unitPathL, unitPathR;

    unitPathL=masterNetwork.getBoundarySection(activeUnit,false);
    unitPathR=masterNetwork.getBoundarySection(activeUnit,true);

    //for now pick longest if not longer than the basic ring size that can be made
    if(unitPathL.size()>=basicMaxSize) unitPath=unitPathR;
    else if(unitPathR.size()>=basicMaxSize) unitPath=unitPathL;
    else if(unitPathL.size()>unitPathR.size()) unitPath=unitPathL;
    else unitPath=unitPathR;

    return unitPath;
}

template <typename CrdT, typename NetT>
void Simulation<CrdT,NetT>::addBasicRing(vector<int> unitPath) {
    //calculate energy of adding basic rings to network, and select by monte carlo method

    vector<int> trialSizes;
    vector<double> trialEnergies;
    trialSizes.clear();
    trialEnergies.clear();
    for(int i=0, j=basicMinSize; i<nBasicRingSizes; ++i, ++j){
        if(j>unitPath.size()){
            masterNetwork.trialRing(j,unitPath,potentialModel);
            trialSizes.push_back(j);
            trialEnergies.push_back(masterNetwork.getEnergy());
        }
    }
    int acceptedRing=monteCarlo.metropolis(trialEnergies);
    if(trialEnergies[acceptedRing]>energyCutoff) killGrowth=true;
    int acceptedSize=trialSizes[acceptedRing];
    masterNetwork.acceptRing(acceptedSize,unitPath,potentialModel);
}
//template <typename CrdT, typename NetT>
//Simulation<CrdT,NetT>::
//template <typename CrdT, typename NetT>
//Simulation<CrdT,NetT>::
//template <typename CrdT, typename NetT>
//Simulation<CrdT,NetT>::
//##### ANALYSE #####
template <typename CrdT, typename NetT>
void Simulation<CrdT,NetT>::analyseNetwork(Logfile &logfile) {
    //analyse network properties
    logfile.log("Analysing network","","",0,false);

    //check geometry
    masterNetwork.checkOverlap();
    logfile.log("Network checked for unit overlap","","",1,false);

    //ring statistics
    masterNetwork.calculateRingStatistics();
    logfile.log("Ring statistics calculated","","",1,false);

    logfile.log("Analysis complete","","",0,true);
}

//##### WRITE #####
template <typename CrdT, typename NetT>
void Simulation<CrdT,NetT>::writeNetwork(Logfile &logfile) {
    //write network to file
    masterNetwork.write(prefixOut,logfile);
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
//template <typename CrdT, typename NetT>
//Simulation<CrdT,NetT>::
//template <typename CrdT, typename NetT>
//Simulation<CrdT,NetT>::
