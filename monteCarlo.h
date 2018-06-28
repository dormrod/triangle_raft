#ifndef MX2_MONTECARLO_H
#define MX2_MONTECARLO_H

#include <iostream>
#include <random>
#include "logfile.h"

using namespace std;

class MonteCarlo {
    //evaluate metropolis condition
private:
    mt19937 mtGen; //mersenne twister generator
    uniform_real_distribution<double> rand01; //uniform distribution
    double rTemperature; //reciprocal temperature

public:
    //constructors
    MonteCarlo();
    MonteCarlo(int seed, double temperature, Logfile &logfile);

};


#endif //MX2_MONTECARLO_H
