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
void Simulation<CrdT,NetT>::setPM(double kMX, double r0MX, double kXX, double a0XX, double kMM, double a0MM, double kLJ, double r0LJ, double kC, double r0C, Logfile &logfile) {
    //set potential model as M-X, X-X, M-M k and r0 values for harmonics
    potentialModel.resize(10);
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
    //geometry specific parameters
    potentialModel[8]=kC;
    potentialModel[9]=r0C;
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
    if(dimensionality==2) masterNetwork=NetT(prefixIn,logfile);
    else if(dimensionality==3) masterNetwork=NetT(prefixIn,logfile,potentialModel[9]);
    masterNetwork.setGO(goMaxIterations,goLineSeachInc,goConvergence,goLocalExtent);
    if(globalPreGO) masterNetwork.geometryOptimiseGlobal(potentialModel);
}

//##### GROW #####
template <typename CrdT, typename NetT>
void Simulation<CrdT,NetT>::growNetwork(Logfile &logfile) {
    //grow network using monte carlo process
    logfile.log("Monte Carlo process","","",0,false);
    logfile.log("Network growth begun after","","sec",1,false);

    //initialise variables
    killGrowth=false;
    energyCutoff=100.0*potentialModel[0];
    goMonitoring=col_vector<int>(3); //number of minimisations, total number of iterations, number of times iteration limit reached
    mcMonitoring=col_vector<double>(nBasicRingSizes);
    int nRings=masterNetwork.getNRings();

    //main loop
    if(nRings<nTargetRings){//only if network needs growing
        do{
            int activeUnit = selectActiveUnit();
            vector<int> unitPath = selectUnitPath(activeUnit);
            addBasicRing(unitPath);
            ++nRings;
            if(nRings%100==0){
                cout<<nRings<<endl;
                logfile.log(to_string(nRings)+" rings, time elapsed: ","","sec",2,false);
            }
            cout<<nRings<<endl;
            if(killGrowth){
                masterNetwork.writeNetwork(prefixOut,logfile);
                logfile.errorlog("Growth prematurely killed due to excessive energy","critical");
                break;
            }
        }while(nRings<nTargetRings);
    }
    logfile.log("All rings built, time elapsed: ","","sec",2,false);
    logfile.log("Network growth complete","","",1,false);

    //monitoring results
    logfile.log("Monte Carlo monitoring","","",1,false);
    logfile.log("Ring proposal probabilities","","",2,false);
    vector<int> ringSizes(nBasicRingSizes);
    for(int i=0, j=basicMinSize; i<nBasicRingSizes; ++i, ++j) ringSizes[i]=j;
    mcMonitoring/=mcMonitoring.sum();
    bool warning=false;
    for(int i=0; i<nBasicRingSizes; ++i) if(mcMonitoring[i]<0.75/nBasicRingSizes) warning=true;
    logfile.log(ringSizes,3,false);
    logfile.log(mcMonitoring,3,false);
    if(warning) logfile.log("Warning: large discrepancy in move proposal probability","","",2,false);
    else logfile.log("Ring proposal probabilities within expected range","","",2,false);
    logfile.log("Geometry optimisation monitoring","","",2,false);
    logfile.log("Average iterations: ",double(goMonitoring[1])/double(goMonitoring[0]),"",3,false);
    logfile.log("Number of times iteration limit reached: ",goMonitoring[2],"",3,false);
    warning=false;
    if(double(goMonitoring[2])/double(goMonitoring[0])>0.05) warning=true;
    if(warning) logfile.log("Warning: iteration limit reached frequently","","",2,false);
    else logfile.log("Iteration limit satisfactory","","",2,false);
    logfile.log("Monitoring analysis complete","","",1,false);

    //global geometry optimisation
    if(globalPostGO){
        logfile.log("Performing global geometry optimisation","","",1,false);
        masterNetwork.geometryOptimiseGlobal(potentialModel);
    }

    logfile.log("Monte Carlo process complete","","",0,true);
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
            //trial ring
            masterNetwork.trialRing(j,unitPath,potentialModel);
            trialSizes.push_back(j);
            trialEnergies.push_back(masterNetwork.getEnergy());
            //monitoring
            mcMonitoring[i]+=1.0;
            int iterations=masterNetwork.getIterations();
            ++goMonitoring[0];
            goMonitoring[1]+=iterations;
            if(iterations==goMaxIterations) ++goMonitoring[2];
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
