#include <iostream>
#include "logfile.h"
#include "simulation.h"

using namespace std;

int main(){

    //Setup logfile and make header
    Logfile logfile("mx2");
    logfile.log("Grows MX2 Aperiodic Networks","","none",0,false);
    logfile.log("Written By: David OM, Wilson Group, 2018","","",0,true);
    logfile.log("Email bugs with snapshots and output files to: ","","",0,false);
    logfile.log("david.ormrodmorley@chem.ox.ac.uk","","",0,true);

    //Open and read main input file
    ifstream inputFile("mx2.inpt", ios::in);
    if(!inputFile.good()) logfile.errorlog("Cannot find input file","critical"); //exit if cannot find input file

    //IO
    string inputPrefix,outputPrefix;
//    int rstFrequency;
    readFileSkipLines(inputFile); //skip header
    readFileValue(inputFile,inputPrefix); //for reading in files
    readFileValue(inputFile,outputPrefix); //for writing to files
//    readFileValue(inputFile,rstFrequency); //restart write out frequency
    //Network properties
    int nTotalRings;
    string geometry;
    int minBasicRingSize, maxBasicRingSize;
    int dimensionality;
    vector<int> basicRingSizeLimits;
    readFileSkipLines(inputFile,2); //skip
    readFileValue(inputFile,nTotalRings); //total rings in generated sample
    readFileRowVector(inputFile,basicRingSizeLimits,2);
    readFileValue(inputFile,geometry); //geometry code for sample
    minBasicRingSize=basicRingSizeLimits[0];
    maxBasicRingSize=basicRingSizeLimits[1];
    if(geometry.substr(0,2)=="2D") dimensionality=2;
    else if(geometry.substr(0,2)=="3D") dimensionality=3;
    //MC
    int randomSeed;
    double temperature;
    readFileSkipLines(inputFile,2); //skip
    readFileValue(inputFile,randomSeed); //for random number generator
    readFileValue(inputFile,temperature); //for metropolis criteria
    //Potential
    double kMX, r0MX, kXX, a0XX, kMM, a0MM, kLJ, r0LJ, kC, r0C;
    vector<double> potential;
    readFileSkipLines(inputFile,2); //skip
    readFileRowVector(inputFile,potential,2); //M-X
    kMX=potential[0];
    r0MX=potential[1];
    readFileRowVector(inputFile,potential,2); //M-X
    kXX=potential[0];
    a0XX=potential[1];
    readFileRowVector(inputFile,potential,2); //M-X
    kMM=potential[0];
    a0MM=potential[1];
    readFileRowVector(inputFile,potential,2); //L-J
    kLJ=potential[0];
    r0LJ=potential[1];
    readFileRowVector(inputFile,potential,2); //Constraints
    kC=potential[0];
    r0C=potential[1];
    //Minimisation
    bool preOpt, postOpt;
    vector<bool> globalOpt;
    int maxIt, localSize;
    double lsInc, convTest;
    readFileSkipLines(inputFile,2); //skip
    readFileRowVector(inputFile,globalOpt,2); //perform global optimisation before/after simulation
    readFileValue(inputFile,maxIt); //maximum iterations of geometry optimisation
    readFileValue(inputFile,lsInc); //line search increment
    readFileValue(inputFile,convTest); //convergence test
    readFileValue(inputFile,localSize); //size of local region
    preOpt=globalOpt[0];
    postOpt=globalOpt[1];
    inputFile.close();

    //set up simulation of correct geometry
    if(geometry=="2DC") {
        Simulation<Cart2D, NetworkCart2D> simulation(logfile);
        simulation.setIO(inputPrefix, outputPrefix, logfile);
        simulation.setNP(nTotalRings, minBasicRingSize, maxBasicRingSize, geometry, logfile);
        simulation.setMC(randomSeed, temperature, logfile);
        simulation.setPM(kMX, r0MX, kXX, a0XX, kMM, a0MM, kLJ, r0LJ, kC, r0C, logfile);
        simulation.setGO(preOpt, postOpt, maxIt, lsInc, convTest, localSize, logfile);

        //run simulation
        simulation.run(logfile);
    }
//    else if(geometry=="3DS") {
//        Simulation<Cart3D, NetworkCart3DS> simulation(logfile);
//        simulation.setIO(inputPrefix, outputPrefix, logfile);
//        simulation.setNP(nTotalRings, minBasicRingSize, maxBasicRingSize, geometry, logfile);
//        simulation.setMC(randomSeed, temperature, logfile);
//        simulation.setPM(kMX, r0MX, kXX, a0XX, kMM, a0MM, kLJ, r0LJ, kC, r0C, logfile);
//        simulation.setGO(preOpt, postOpt, maxIt, lsInc, convTest, localSize, logfile);
//
//        //run simulation
//        simulation.run(logfile);
//    }

    return 0;
}
