//Geometry optimisation functors that take a coordinate type and potential model as template parameters
#ifndef MX2_GEOM_OPT_ALGS_H
#define MX2_GEOM_OPT_ALGS_H

#include <iostream>
#include <limits>
#include "col_vector.h"

using namespace std;

template <typename PotModel>
class SteepestDescent{
private:
    int iterationLimit; //maximum iterations
    double lineInc; //line search increment
    double convCriteria; //convergence criteria

public:
    //constructors
    SteepestDescent();
    SteepestDescent(int maxIt, double lInc, double cc);

    //setters
    void setLineSearchIncrement(double inc);

    //function call
    int operator()(PotModel &model, double &energy, int &iterations, vector<double> &crdsIn);
};

template <typename PotModel>
class SteepestDescentArmijo{
private:
    int iterationLimit; //maximum iterations
    double tau; //line search increment
    double convCriteria; //convergence criteria

public:
    //constructors
    SteepestDescentArmijo();
    SteepestDescentArmijo(int maxIt, double t, double cc);

    //function call
    int operator()(PotModel &model, double &energy, int &iterations, vector<double> &crdsIn);
};

#include "geom_opt_algs.tpp"
#endif //MX2_GEOM_OPT_ALGS_H
