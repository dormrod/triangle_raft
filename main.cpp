#include <iostream>
#include "logfile.h"
#include "netBaseUnits.h"

using namespace std;

int main(){

    //Setup logfile and make header
    Logfile logfile("mx2");
    logfile.log("Grows MX2 Aperiodic Networks","","none",0,false);
    logfile.log("Written By: David OM, Wilson Group, 2018","","",0,true);
    logfile.log("Email bugs with snapshots and output files to: ","","",0,false);
    logfile.log("david.ormrodmorley@chem.ox.ac.uk","","",0,true);

    //Open and read main input file
//    ifstream inputFile("tang.inpt", ios::in);
//    if(!inputFile.good()) logfile.errorlog("Cannot find input file","critical"); //exit if cannot find input file

    return 0;
}
