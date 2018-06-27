#ifndef MX2_NETWORKDERIVED_H
#define MX2_NETWORKDERIVED_H

#include <iostream>
#include <string>
#include "logfile.h"
#include "crd.h"
#include "networkBase.h"

using namespace std;

class NetworkCart2D: public Network<Cart2D> {
    //network class using two dimensional cartesian coordinate
private:

public:
    //constructors
    NetworkCart2D(string prefix, Logfile &logfile); //load network from files

    //write out
    void write(string prefix, Logfile &logfile) override; //write network to files
};


#endif //MX2_NETWORKDERIVED_H
