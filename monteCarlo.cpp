#include "monteCarlo.h"

MonteCarlo::MonteCarlo() {
    //default constructor - initialise with seed 0
    mtGen.seed(0); //mersenne twister generator
    rand01=uniform_real_distribution<double>(0.0,1.0); //uniform distribution between 0->1
    rTemperature=1.0;
}

MonteCarlo::MonteCarlo(int seed, double temperature, Logfile &logfile) {
    //constructor - intitialise with given seed
    mtGen.seed(seed); //mersenne twister generator
    rand01=uniform_real_distribution<double>(0.0,1.0); //uniform distribution between 0->1
    rTemperature=1.0/temperature;
    logfile.log("Initialised: ","Monte Carlo, with seed "+to_string(seed),"",1,false);
}

int MonteCarlo::metropolis(vector<double> energies) {
    //apply metropolis condition to relative energies
    int n = energies.size();
    double lowestEnergy = numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) if (energies[i] < lowestEnergy) lowestEnergy = energies[i];
    vector<double> deltaE = energies;
    for (int i = 0; i < n; ++i) deltaE[i] -= lowestEnergy;
    vector<double> probabilities;
    double p, pTot = 0.0;
    for (int j = 0; j < n; ++j) {
        p = exp(-deltaE[j] * rTemperature);
        probabilities.push_back(p);
        pTot += p;
    }
    for (int j = 0; j < n; ++j) probabilities[j] /= pTot;
    double r01 = rand01(mtGen);
    int acceptedNetwork;
    pTot = 0.0;
    for (int j = 0; j < n; ++j) {
        pTot += probabilities[j];
        if (r01 < pTot) {
            acceptedNetwork = j;
            break;
        }
    }
    return acceptedNetwork;
}