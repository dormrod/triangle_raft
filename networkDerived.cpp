#include "networkDerived.h"

//##### NETWORK CART 2D #####
NetworkCart2D::NetworkCart2D() {
    //default constructor
}

NetworkCart2D::NetworkCart2D(string prefix, Logfile &logfile):Network<Cart2D>() {
    //load network from files, except coordinate data

    //set up file names
    string atomFilename=prefix+"_atoms.out";
    string unitFilename=prefix+"_units.out";
    string ringFilename=prefix+"_rings.out";
    string cnxFilename=prefix+"_connections.out";

    //check files data
    logfile.log("Loading network","","",0,false);
    ifstream atomFile(atomFilename, ios::in);
    ifstream unitFile(unitFilename, ios::in);
    ifstream ringFile(ringFilename, ios::in);
    ifstream cnxFile(cnxFilename, ios::in);
    if(!atomFile.good()) logfile.errorlog("Cannot find atom input file","critical");
    if(!unitFile.good()) logfile.errorlog("Cannot find unit input file","critical");
    if(!ringFile.good()) logfile.errorlog("Cannot find ring input file","critical");
    if(!cnxFile.good()) logfile.errorlog("Cannot find connection input file","critical");

    //read data
    vector< vector<int> > dataI;
    vector< vector<double> > dataD;

    //set up all atoms
    readFileAll(atomFile,dataD);
    for(int i=0; i<dataD.size();++i){
        int elem=int(dataD[i][0]);
        int cnd=int(dataD[i][1]);
        Cart2D crd = Cart2D(dataD[i][2], dataD[i][3]);
        Atom<Cart2D> atom(i,elem,cnd);
        atom.setCoordinate(crd);
        addAtom(atom);
    }

    //inital set up of units
    readFileAll(unitFile,dataI);
    for(int i=0; i<dataI.size();++i){//set up
        Unit triangle(i,3,3,3);
        addUnit(triangle);
    }
    for(int i=0; i<dataI.size();++i){//add m atom and unit-atom connections
        int m=dataI[i][0];
        units[i].setAtomM(m);
        for(int j=0; j<3; ++j){
            addUnitAtomXCnx(i,dataI[i][j+1]);
        }
    }

    //initial set up of rings
    readFileAll(ringFile,dataI);
    for(int i=0; i<dataI.size();++i){//set up
        Ring ring(i,dataI[i].size(),dataI[i].size());
        addRing(ring);
    }
    for(int i=0; i<dataI.size();++i){//add unit-ring connections
        for(int j=0; j<dataI[i].size(); ++j){
            addUnitRingCnx(dataI[i][j],i);
        }
    }

    //add unit-unit connections
    int n;
    vector<int> pair;
    readFileValue(cnxFile,n);
    for(int i=0; i<n; ++i){
        readFileRowVector(cnxFile,pair,2);
        addUnitUnitCnx(pair[0],pair[1]);
    }
    //add ring-ring connections
    readFileValue(cnxFile,n);
    for(int i=0; i<n; ++i){
        readFileRowVector(cnxFile,pair,2);
        addRingRingCnx(pair[0],pair[1]);
    }

    //close files
    atomFile.close();
    unitFile.close();
    ringFile.close();
    cnxFile.close();

    logfile.log("Atoms read from: ", atomFilename, "", 1, false);
    logfile.log("Units read from: ", unitFilename, "", 1, false);
    logfile.log("Rings read from: ", ringFilename, "", 1, false);
    logfile.log("Additional connections read from: ", cnxFilename, "", 1, false);
    logfile.log("Load complete","","",0,true);
}

void NetworkCart2D::write(string prefix, Logfile &logfile) {
    //write network information to output files

    //set up file names
    string atomFilename=prefix+"_atoms.out";
    string unitFilename=prefix+"_units.out";
    string ringFilename=prefix+"_rings.out";
    string cnxFilename=prefix+"_connections.out";

    //write out all atom, unit, ring and connection data
    logfile.log("Writing network","","",0,false);
    ofstream atomFile(atomFilename, ios::in|ios::trunc);
    ofstream unitFile(unitFilename, ios::in|ios::trunc);
    ofstream ringFile(ringFilename, ios::in|ios::trunc);
    ofstream cnxFile(cnxFilename, ios::in|ios::trunc);

    //write atom element, coordination and x-y coordinate
    atomFile << fixed << showpoint << setprecision(6);
    for(int i=0; i<nAtoms; ++i){
        atomFile<<atoms[i].element<<"  "<<atoms[i].coordination<<"  "<<atoms[i].coordinate.x<<"  "<<atoms[i].coordinate.y<<endl;
    }
    logfile.log("Atoms written to: ", atomFilename, "", 1, false);

    //write unit atom ids
    for(int i=0; i<nUnits; ++i){
        writeFileValue(unitFile,units[i].atomM,false);
        writeFileArray(unitFile,units[i].atomsX.ids,3,true);
    }
    logfile.log("Units written to: ", unitFilename, "", 1, false);

    //write ring unit geometrical unit ids
    for(int i=0; i<nRings; ++i){
        writeFileArray(ringFile,rings[i].units.ids,rings[i].units.n,false);
    }
    logfile.log("Rings written to: ", ringFilename, "", 1, false);

    //write additional connections
    int p0, p1;
    vector<int> pair(2);
    vector< vector<int> > unitUnitCnxs;
    unitUnitCnxs.clear();
    for(int i=0; i<nUnits; ++i){
        p0=units[i].id;
        for(int j=0; j<units[i].units.n; ++j){
            p1=units[units[i].units.ids[j]].id;
            if(p0<p1){
                pair[0]=p0;
                pair[1]=p1;
                unitUnitCnxs.push_back(pair);
            }
        }
    }
    writeFileValue(cnxFile,unitUnitCnxs.size(),true);
    for(int i=0; i<unitUnitCnxs.size(); ++i) writeFileVector(cnxFile,unitUnitCnxs[i]);
    vector< vector<int> > ringRingCnxs;
    ringRingCnxs.clear();
    for(int i=0; i<nRings; ++i){
        p0=rings[i].id;
        for(int j=0; j<rings[i].rings.n; ++j){
            p1=rings[rings[i].rings.ids[j]].id;
            if(p0<p1){
                pair[0]=p0;
                pair[1]=p1;
                ringRingCnxs.push_back(pair);
            }
        }
    }
    writeFileValue(cnxFile,ringRingCnxs.size(),true);
    for(int i=0; i<ringRingCnxs.size(); ++i) writeFileVector(cnxFile,ringRingCnxs[i]);

    logfile.log("Additional connections written to: ", cnxFilename, "", 1, false);

    atomFile.close();
    unitFile.close();
    ringFile.close();
    cnxFile.close();
    logfile.log("Write complete","","",0,true);
}

void NetworkCart2D::setGO(int it, double ls, double conv) {
    //set up optimiser with geometry optimisation parameters
    optimiser=SteepestDescent<HC2>(it,ls,conv);
}

vector<double> NetworkCart2D::getCrds() {
    //get all atom coordinates collapsed onto 1d
    vector<double> crds;
    crds.clear();
    for(int i=0; i<nAtoms; ++i){
        crds.push_back(atoms[i].coordinate.x);
        crds.push_back(atoms[i].coordinate.y);
    }
    return crds;
}

void NetworkCart2D::setCrds(vector<double> &crds) {
    //set all atom coordinates
    for(int i=0; i<nAtoms; ++i){
        atoms[i].coordinate.x=crds[2*i];
        atoms[i].coordinate.y=crds[2*i+1];
    }
}

void NetworkCart2D::geometryOptimise(vector<double> &potentialModel) {
    //set up harmonic potential before passing to derived class

    //reset potential information - don't need angles for harmonic potential
    vector<int> bonds, angles, fixedAtoms, interx;
    vector<double>  bondK, bondR0, angleK, angleR0, crds;
    bonds.clear();
    angles.clear();
    fixedAtoms.clear();
    interx.clear();
    bondK.clear();
    bondR0.clear();
    angleK.clear();
    angleR0.clear();

    //get all atom coordinates
    crds=getCrds();

    //loop over triangle units to get M-X, X-X bonds
    int mId0, xId0, xId1, xId2;
    for(int i=0; i<nUnits; ++i){
        mId0=units[i].atomM;
        xId0=units[i].atomsX.ids[0];
        xId1=units[i].atomsX.ids[1];
        xId2=units[i].atomsX.ids[2];
        //M-X
        bonds.push_back(mId0);
        bonds.push_back(xId0);
        bondK.push_back(potentialModel[0]);
        bondR0.push_back(potentialModel[1]);
        bonds.push_back(mId0);
        bonds.push_back(xId1);
        bondK.push_back(potentialModel[0]);
        bondR0.push_back(potentialModel[1]);
        bonds.push_back(mId0);
        bonds.push_back(xId2);
        bondK.push_back(potentialModel[0]);
        bondR0.push_back(potentialModel[1]);
        //X-X
        bonds.push_back(xId0);
        bonds.push_back(xId1);
        bondK.push_back(potentialModel[2]);
        bondR0.push_back(potentialModel[3]);
        bonds.push_back(xId0);
        bonds.push_back(xId2);
        bondK.push_back(potentialModel[2]);
        bondR0.push_back(potentialModel[3]);
        bonds.push_back(xId1);
        bonds.push_back(xId2);
        bondK.push_back(potentialModel[2]);
        bondR0.push_back(potentialModel[3]);
    }

    //loop over neighbour triangle units to get M-M
    int mId1;
    for(int i=0; i<nUnits; ++i){
        mId0=units[i].atomM;
        for(int j=0; j<units[i].units.n; ++j){
            mId1=units[units[i].units.ids[j]].atomM;
            if(mId0<mId1){//prevent double counting as reciprocal connections
                bonds.push_back(mId0);
                bonds.push_back(mId1);
                bondK.push_back(potentialModel[4]);
                bondR0.push_back(potentialModel[5]);
            }
        }
    }

    //set up model and optimise
    HC2 harmonicPotential(bonds, angles, bondK, bondR0, angleK, angleR0, fixedAtoms, interx);
    optimiser(harmonicPotential, energy, optIterations, crds);

    //update coordinates
    setCrds(crds);
}
