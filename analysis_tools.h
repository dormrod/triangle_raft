#ifndef MX2_ANALYSIS_TOOLS_H
#define MX2_ANALYSIS_TOOLS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include "col_vector.h"
#include "easyIO.h"

using namespace std;

struct DiscreteDistribution{
    //holds discrete probability distribution and statistics

    //distribution variables
    int n; //number of values
    col_vector<int> x; //values
    col_vector<double> p, p_raw; //(normalised) probabilities

    //statistics
    double mean;

    //constructors, destructors
    DiscreteDistribution(); //default
    DiscreteDistribution(vector<int> values);

    //getters
    vector<int> getValues();
    vector<double> getProbabilities();
    double getProbability(int xValue);

};

struct ContinuousDistribution{
    //holds continuous probability distribution and statistics

    //distribution variables
    int n; //number of values
    col_vector<double> x; //values

    //statistics
    double mean;

    //constructors, destructors
    ContinuousDistribution(); //default
    ContinuousDistribution(vector<double> values);

};


//standalone functions
vector<double> leastSquaresLinearRegression(vector<double> x, vector<double> y);

#endif //MX2_ANALYSIS_TOOLS_H
