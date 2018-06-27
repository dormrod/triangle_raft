#include "networkDerived.h"

//##### NETWORK CART 2D #####

NetworkCart2D::NetworkCart2D(string prefix, Logfile &logfile):Network<Cart2D>() {
    //load network from files, except coordinate data

    //set up file names
    string atomFilename=prefix+"_atoms.out";
    string unitFilename=prefix+"_units.out";
    string ringFilename=prefix+"_rings.out";

    //read in data
    logfile.log("Loading network","","",0,false);
    ifstream atomFile(atomFilename, ios::in);
    ifstream unitFile(unitFilename, ios::in);
    ifstream ringFile(ringFilename, ios::in);
    //atom coordinates
    if(!atomFile.good()) logfile.errorlog("Cannot find atom input file","critical");
    else {
        vector<vector<double> > atomsIn;
        readFileAll(atomFile, atomsIn);
        nAtoms = atomsIn.size();
        atoms.resize(nAtoms);
        int elem; //element
        Cart2D crd; //x,y coordinate
        for (int i = 0; i < nAtoms; ++i) {
            elem = int(atomsIn[i][0]);
            crd = Cart2D(atomsIn[i][1], atomsIn[i][2]);
            Atom<Cart2D> atom(i, elem);
            atom.setCoordinate(crd);
            atoms[i] = atom;
        }
        logfile.log("Atoms read from: ", atomFilename, "", 1, false);
    }
    //triangular geometrical units
    if(!unitFile.good()) logfile.errorlog("Cannot find unit input file","critical");
    else{
        vector< vector<int> > geomUnitsIn;
        readFileAll(unitFile,geomUnitsIn);
        nGeomUnits=geomUnitsIn.size();
        geomUnits.resize(nGeomUnits);
        int m; //metal atom as last in list
        vector<int> x(3); //x atoms as first in list
        for(int i=0; i<nGeomUnits; ++i){
            m=geomUnitsIn[i][3];
            for(int j=0; j<3; ++j) x[j]=geomUnitsIn[i][j];
            GeometricalUnit<Cart2D> triangle(i,3);
            triangle.setAtomM(m);
            triangle.setAtomsX(x);
            geomUnits[i]=triangle;
        }
        logfile.log("Geometrical units read from: ",unitFilename,"",1,false);
    }
    //rings of units
    if(!ringFile.good()) logfile.errorlog("Cannot find ring input file","critical");
    else{
        vector< vector<int> > ringsIn;
        readFileAll(ringFile,ringsIn);
        nRings=ringsIn.size();
        rings.resize(nRings);
        for(int i=0; i<nRings; ++i){
            Ring<Cart2D> ring(i,ringsIn[i].size());
            ring.setUnits(ringsIn[i]);
            rings[i]=ring;
        }
        logfile.log("Rings read from: ",ringFilename,"",1,false);
    }
    atomFile.close();
    unitFile.close();
    ringFile.close();
    logfile.log("Load complete","","",0,true);
}

void NetworkCart2D::write(string prefix, Logfile &logfile) {
    //write network information to output files

    //set up file names
    string atomFilename=prefix+"_atoms.out";
    string unitFilename=prefix+"_units.out";
    string ringFilename=prefix+"_rings.out";

    //write out all atom, geometrical unit and ring data
    logfile.log("Writing network","","",0,false);
    ofstream atomFile(atomFilename, ios::in|ios::trunc);
    ofstream unitFile(unitFilename, ios::in|ios::trunc);
    ofstream ringFile(ringFilename, ios::in|ios::trunc);

    //write atom element and x-y coordinate
    atomFile << fixed << showpoint << setprecision(6);
    for(int i=0; i<nAtoms; ++i){
        atomFile<<atoms[i].element<<"  "<<atoms[i].coordinate.x<<"  "<<atoms[i].coordinate.y<<endl;
    }
    logfile.log("Atoms written to: ", atomFilename, "", 1, false);

    //write geometrical unit atom ids
    for(int i=0; i<nGeomUnits; ++i){
        writeFileArray(unitFile,geomUnits[i].xAtoms,3,false);
        writeFileValue(unitFile,geomUnits[i].mAtom);
    }
    logfile.log("Units written to: ", unitFilename, "", 1, false);

    //write ring unit geometrical unit ids
    for(int i=0; i<nRings; ++i){
        writeFileArray(ringFile,rings[i].units,rings[i].n,false);
    }
    logfile.log("Rings written to: ", ringFilename, "", 1, false);

    atomFile.close();
    unitFile.close();
    ringFile.close();
    logfile.log("Write complete","","",0,true);
}