#include "geom_opt_algs.h"

//##### STEEPEST DESCENT #####
template <typename PotModel>
SteepestDescent<PotModel>::SteepestDescent() {
    //default constructor
    iterationLimit=0;
}

template <typename PotModel>
SteepestDescent<PotModel>::SteepestDescent(int maxIt, double lInc, double cc) {
    //set steepest descent parameters
    iterationLimit=maxIt;
    lineInc=lInc;
    convCriteria=cc;
}

template <typename PotModel>
void SteepestDescent<PotModel>::setLineSearchIncrement(double inc) {
    //alter line search increment
    lineInc=inc;
}

template <typename PotModel>
int SteepestDescent<PotModel>::operator()(PotModel &model, double &energy, int &iterations, vector<double> &crdsIn) {
    //steepest descent algorithm
    col_vector<double> crds=crdsIn;
    col_vector<double> force(crds.n);

    //intialise steepest descent parameters
    energy=0.0;
    iterations=0;
    double previousEnergy=numeric_limits<double>::infinity();
    double deltaE; //difference between current and previous energy
    col_vector<double> crdInc; //increment in coordinates for line search

    //evaluate force and check non-zero before commencing main loop
    model.calculateForce(crds,force);
    if(force.asum()<1e-6) return 1;

    //steepest descent algorithm
    for(int i=0; i<iterationLimit; ++i){
        //line search
        double e0, e1;
        crdInc=force*lineInc;
        model.calculateEnergy(crds,e0);
        for(;;){
            crds+=crdInc;
            model.calculateEnergy(crds,e1);
            if(e1>e0){//passed through minimum so backtrack
                crds-=crdInc;
                energy=e0;
                ++iterations;
                break;
            }
            else e0=e1;
        }

        //check energy convergence
        deltaE=fabs(energy-previousEnergy);
        if(deltaE<convCriteria) break;
        else previousEnergy=energy;

        //recalculate forces
        model.calculateForce(crds,force);
//        cout<<i<<" "<<energy<<" "<<deltaE<<endl;
    }

//    cout<<"iterations "<<iterations<<" energy "<<energy<<endl;

    //update coordinates
    for(int i=0; i<crds.n; ++i) crdsIn[i]=crds[i];
    return 0;
}