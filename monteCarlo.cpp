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