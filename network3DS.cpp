#include "network3DS.h"

//##### NETWORK CART 3D S #####
NetworkCart3DS::NetworkCart3DS() {
    //default constructor
}

NetworkCart3DS::NetworkCart3DS(string prefix, Logfile &logfile, double additionalParams):Network<Cart3D>() {
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
        Cart3D crd = Cart3D(dataD[i][2], dataD[i][3], additionalParams); //set z coordinate as initial sphere radius
        Atom<Cart3D> atom(i,elem,cnd);
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

    //calculate boundary
    calculateBoundary();

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

void NetworkCart3DS::writeNetwork(string prefix, Logfile &logfile) {
    //write network information to output files

    //set up file names
    string atomFilename=prefix+"_atoms.out";
    string unitFilename=prefix+"_units.out";
    string ringFilename=prefix+"_rings.out";
    string cnxFilename=prefix+"_connections.out";
    string visFilename=prefix+"_vis.out";
    string xyzFilename=prefix+".xyz";

    //write out all atom, unit, ring and connection data
    logfile.log("Writing network","","",0,false);
    ofstream atomFile(atomFilename, ios::in|ios::trunc);
    ofstream unitFile(unitFilename, ios::in|ios::trunc);
    ofstream ringFile(ringFilename, ios::in|ios::trunc);
    ofstream cnxFile(cnxFilename, ios::in|ios::trunc);
    ofstream visFile(visFilename, ios::in|ios::trunc);

    //write atom element, coordination and x,y,z coordinate
    atomFile << fixed << showpoint << setprecision(6);
    for(int i=0; i<nAtoms; ++i){
        atomFile<<atoms[i].element<<"  "<<atoms[i].coordination<<"  "<<atoms[i].coordinate.x<<"  "<<atoms[i].coordinate.y<<" "<<atoms[i].coordinate.z<<endl;
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
        writeFileArray(ringFile,rings[i].units.ids,rings[i].units.n,true);
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

    //write ring colour codes
    for(int i=0; i<nRings; ++i) writeFileVector(visFile,ringColours[i]);
    logfile.log("Visualisation helper file written to: ", visFilename, "", 1, false);

    //write x,y,z file
    ofstream xyzFile(xyzFilename, ios::in|ios::trunc);
    atomFile << fixed << showpoint << setprecision(6);
    writeFileValue(xyzFile,nAtoms,true);
    writeFileValue(xyzFile," ",true);
    for(int i=0; i<nAtoms; ++i){
        xyzFile<<atoms[i].element<<"  "<<atoms[i].coordinate.x<<"  "<<atoms[i].coordinate.y<<" "<<atoms[i].coordinate.z<<endl;
    }
    logfile.log("xyz file written to: ", xyzFilename, "", 1, false);

    atomFile.close();
    unitFile.close();
    ringFile.close();
    cnxFile.close();
    visFile.close();
    xyzFile.close();
    logfile.log("Write complete","","",0,true);
}

void NetworkCart3DS::writeNetworkSpecial(string prefix, Logfile &logfile) {
    //placeholder
    return;
}

vector<double> NetworkCart3DS::getCrds() {
    //get all atom coordinates collapsed onto 1d
    vector<double> crds;
    crds.clear();
    for(int i=0; i<nAtoms; ++i){
        crds.push_back(atoms[i].coordinate.x);
        crds.push_back(atoms[i].coordinate.y);
        crds.push_back(atoms[i].coordinate.z);
    }
    return crds;
}

vector<double> NetworkCart3DS::getCrds(map<int, int> &globalAtomMap, int n) {
    //get local atom coordinates collapsed onto 1d
    vector<double> crds;
    crds.clear();
    int id;
    for(int i=0; i<n; ++i){
        id=globalAtomMap.at(i);
        crds.push_back(atoms[id].coordinate.x);
        crds.push_back(atoms[id].coordinate.y);
        crds.push_back(atoms[id].coordinate.z);
    }
    return crds;
}

void NetworkCart3DS::setCrds(vector<double> &crds) {
    //set all atom coordinates
    for(int i=0; i<nAtoms; ++i){
        atoms[i].coordinate.x=crds[3*i];
        atoms[i].coordinate.y=crds[3*i+1];
        atoms[i].coordinate.z=crds[3*i+2];
    }
}

void NetworkCart3DS::setCrds(map<int, int> &globalAtomMap, vector<double> &crds) {
    //set local atom coordinates
    int n=crds.size()/3;
    int id;
    for(int i=0; i<n; ++i){
        id=globalAtomMap.at(i);
        atoms[id].coordinate.x=crds[3*i];
        atoms[id].coordinate.y=crds[3*i+1];
        atoms[id].coordinate.z=crds[3*i+2];
    }
}

void NetworkCart3DS::setGO(int it, double ls, double conv, int loc) {
    //set up optimiser with geometry optimisation parameters
    localExtent=loc;
    defLineInc=ls;
    optimiser=SteepestDescentArmijo<HLJC3S>(it,ls,conv);
}

void NetworkCart3DS::geometryOptimiseGlobal(vector<double> &potentialModel) {
    //set up harmonic potential before passing to derived class

    //reset potential information
    vector<int> bonds, angles, repulsions, fixedAtoms, interx, constrainedAtoms;
    vector<double>  bondK, bondR0, repK, repR0, conK, conR0, crds;
    bonds.clear();
    angles.clear();
    repulsions.clear();
    fixedAtoms.clear();
    interx.clear();
    constrainedAtoms.clear();
    bondK.clear();
    bondR0.clear();
    repK.clear();
    repR0.clear();

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

    //loop over all m atoms and add M-M LJ interactions
    int mId1;
    for(int i=0; i<nUnits-1;++i){
        mId0=units[i].atomM;
        for(int j=i+1; j<nUnits; ++j){
            mId1=units[j].atomM;
            repulsions.push_back(mId0);
            repulsions.push_back(mId1);
            repK.push_back(potentialModel[6]);
            repR0.push_back(potentialModel[7]);
        }
    }

    //add spherical constraints
    for(int i=0; i<nAtoms; ++i){
        constrainedAtoms.push_back(i);
        conK.push_back(potentialModel[8]);
        conR0.push_back(potentialModel[9]);
    }

    //set up model and optimise
    HLJC3S potential(bonds,angles,repulsions,fixedAtoms,interx,constrainedAtoms,bondK,bondR0,repK,repR0,conK,conR0);
    optimiser(potential, energy, optIterations, crds);

    //update coordinates
    setCrds(crds);
}

void NetworkCart3DS::geometryOptimiseLocal(vector<double> &potentialModel) {
    //geometry optimise atoms only in local region

    //reset potential information
    vector<int> bonds, angles, repulsions, fixedAtoms, interx, constrainedAtoms;
    vector<double>  bondK, bondR0, repK, repR0, conK, conR0, crds;
    bonds.clear();
    angles.clear();
    repulsions.clear();
    fixedAtoms.clear();
    interx.clear();
    constrainedAtoms.clear();
    bondK.clear();
    bondR0.clear();
    repK.clear();
    repR0.clear();

    //get local region
    findLocalRegion(rings.rbegin()[0].id,localExtent);

    //get local atom coordinates
    crds=getCrds(globalAtomMap,nLocalAtoms);

    //loop over triangle units to get M-X, X-X bonds
    int mId0, xId0, xId1, xId2;
    for(int i=0; i<flexLocalUnits.size(); ++i){
        mId0=localAtomMap.at(units[flexLocalUnits[i]].atomM);
        xId0=localAtomMap.at(units[flexLocalUnits[i]].atomsX.ids[0]);
        xId1=localAtomMap.at(units[flexLocalUnits[i]].atomsX.ids[1]);
        xId2=localAtomMap.at(units[flexLocalUnits[i]].atomsX.ids[2]);
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
    for(int i=0; i<fixedLocalUnits.size(); ++i){
        mId0=localAtomMap.at(units[fixedLocalUnits[i]].atomM);
        xId0=localAtomMap.at(units[fixedLocalUnits[i]].atomsX.ids[0]);
        xId1=localAtomMap.at(units[fixedLocalUnits[i]].atomsX.ids[1]);
        xId2=localAtomMap.at(units[fixedLocalUnits[i]].atomsX.ids[2]);
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

    //loop over m atoms and add M-M LJ interactions
    int mId1;
    for(int i=0; i<flexLocalUnits.size()-1;++i){
        mId0=localAtomMap.at(units[flexLocalUnits[i]].atomM);
        for(int j=i+1; j<flexLocalUnits.size(); ++j){
            mId1=localAtomMap.at(units[flexLocalUnits[j]].atomM);
            repulsions.push_back(mId0);
            repulsions.push_back(mId1);
            repK.push_back(potentialModel[6]);
            repR0.push_back(potentialModel[7]);
        }
        for(int j=0; j<fixedLocalUnits.size(); ++j){
            mId1=localAtomMap.at(units[fixedLocalUnits[j]].atomM);
            repulsions.push_back(mId0);
            repulsions.push_back(mId1);
            repK.push_back(potentialModel[6]);
            repR0.push_back(potentialModel[7]);
        }
    }

    //add spherical constraints
    for(int i=0; i<nLocalAtoms; ++i){
        constrainedAtoms.push_back(i);
        conK.push_back(potentialModel[8]);
        conR0.push_back(potentialModel[9]);
    }

    //set up model and optimise
    HLJC3S potential(bonds,angles,repulsions,fixedAtoms,interx,constrainedAtoms,bondK,bondR0,repK,repR0,conK,conR0);
    optimiser(potential, energy, optIterations, crds);

    //update coordinates
    setCrds(globalAtomMap,crds);
}

bool NetworkCart3DS::checkGrowth() {
    //check if any boundary atoms in negative z-direction, limits to hemisphere

    int mId;
    double z;
    bool flag=false;
    for(int i=0; i<boundaryUnits.size(); ++i){
        mId=units[boundaryUnits[i]].atomM;
        z=atoms[mId].coordinate.z;
        if(z<0.0) flag=true;
    }
    return flag;
}

bool NetworkCart3DS::checkLocalGrowth(int rId) {
    //placeholder
    return true;
}

int NetworkCart3DS::getActiveUnit(string shape, double size){
//    find active unit within given shape, return -1 if cannot find

    int unitId=-1;
    int atomId; //id of atom with dangling bond

    if(shape=="S") {//projected onto circle in x-y plane
        Cart2D projection;
        double rsq = size * size;
        for (int i = 0; i < boundaryUnits.size(); ++i) {
            atomId = boundaryStatus[i];
            if (atomId >= 0) {
                projection=atoms[atomId].coordinate.xyProjection();
                if (projection.normSq() < rsq) {
                    unitId = boundaryUnits[i];
                    break;
                }
            }
        }
    }
    return unitId;
}

void NetworkCart3DS::buildRing(int ringSize, vector<int> &unitPath, vector<double> &potentialModel){
    //build a ring of a given size to a starting path

    //calculate number of new species
    int nNewTriangles=ringSize-unitPath.size();
    int nNewM=nNewTriangles;
    int nNewX=3*nNewM-1-nNewM;
    int nNewX3=nNewTriangles;
    int nNewX4=nNewX-nNewX3;

    //set up ids for new species
    int ringId=nRings;
    vector<int> mAtomIds(nNewM), x3AtomIds(nNewX3), x4AtomIds(nNewX4), triIds(nNewTriangles);
    for(int i=0; i<nNewTriangles; ++i) triIds[i]=nUnits+i;
    for(int i=0; i<nNewM; ++i) mAtomIds[i]=nAtoms+i;
    for(int i=0; i<nNewX3; ++i) x3AtomIds[i]=mAtomIds.rbegin()[0]+i+1;
    for(int i=0; i<nNewX4; ++i) x4AtomIds[i]=x3AtomIds.rbegin()[0]+i+1;

    //make new atoms
    for(int i=0; i<nNewM; ++i){//m
        Atom<Cart3D> atom(nAtoms,14,3);
        addAtom(atom);
    }
    for(int i=0; i<nNewX3; ++i){//x3
        Atom<Cart3D> atom(nAtoms,8,3);
        addAtom(atom);
    }
    for(int i=0; i<nNewX4; ++i){//x4
        Atom<Cart3D> atom(nAtoms,8,4);
        addAtom(atom);
    }
    //increase coordination of dangling atoms in unit path
    int atomIdL = boundaryStatus[find(boundaryUnits.begin(), boundaryUnits.end(), unitPath[0]) - boundaryUnits.begin()];
    int atomIdR = boundaryStatus[find(boundaryUnits.begin(), boundaryUnits.end(), unitPath.rbegin()[0]) - boundaryUnits.begin()];
    ++atoms[atomIdL].coordination;
    ++atoms[atomIdR].coordination;

    //make new triangles
    for(int i=0; i<nNewTriangles; ++i){
        Unit triangle(nUnits,3,3,3);
        addUnit(triangle);
    }
    //assign m atoms and x3 atoms
    for(int i=0; i<nNewTriangles; ++i){
        units[triIds[i]].setAtomM(mAtomIds[i]);
        addUnitAtomXCnx(triIds[i],x3AtomIds[i]);
    }
    //assign x4 atoms
    for(int i=0; i<nNewTriangles; ++i){
        if(i==0) addUnitAtomXCnx(triIds[i],atomIdL);
        else addUnitAtomXCnx(triIds[i],x4AtomIds[i-1]);
        if(i==nNewX4) addUnitAtomXCnx(triIds[i],atomIdR);
        else addUnitAtomXCnx(triIds[i],x4AtomIds[i]);
    }
    //assign triangles
    addUnitUnitCnx(triIds[0],unitPath[0]);
    for(int i=0; i<nNewTriangles-1; ++i) addUnitUnitCnx(triIds[i],triIds[i+1]);
    addUnitUnitCnx(triIds.rbegin()[0],unitPath.rbegin()[0]);

    //make new ring
    Ring ring(nRings,ringSize,ringSize);
    addRing(ring);
    //assign ring-units
    for(int i=0; i<unitPath.size(); ++i) addUnitRingCnx(unitPath.rbegin()[i],ringId);
    for(int i=0; i<nNewTriangles; ++i) addUnitRingCnx(triIds[i],ringId);
    //assign ring-rings
    vector<int> nbRings;
    nbRings.clear();
    for(int i=0; i<unitPath.size(); ++i){
        for(int j=0; j<units[unitPath[i]].rings.n;++j){
            nbRings.push_back(units[unitPath[i]].rings.ids[j]);
        }
    }
    sort(nbRings.begin(), nbRings.end());
    int prevRing=-1, currRing;
    for(int i=0; i<nbRings.size(); ++i){
        currRing=nbRings[i];
        if(currRing!=prevRing && currRing!=ringId){
            addRingRingCnx(currRing,ringId);
            prevRing=currRing;
        }
    }

    //generate new coordinates by projecting onto x-y plane, then offsetting in z-direction
    //generate vectors
    Cart2D uvPar, uvPer; //unit vectors parallel and perpendicular to L->R
    uvPar=(atoms[atomIdR].coordinate-atoms[atomIdL].coordinate).xyProjection();
    Cart2D vClockwise(uvPar.y,-uvPar.x);
    Cart2D vAntiClockwise(-uvPar.y,uvPar.x);
//    Cart2D dir=atoms[atomIdL].coordinate-atoms[units[unitPath[0]].atomM].coordinate;
    Cart2D dir=((atoms[atomIdL].coordinate+atoms[atomIdR].coordinate)*0.5).xyProjection();
    if(vClockwise*dir>0) uvPer=vClockwise;
    else uvPer=vAntiClockwise;
    double perLen, parLen; //lengths of original perpendicular and parallel vectors
    uvPar.normalise(parLen);
    uvPer.normalise(perLen);

    if (nNewTriangles%2==0){
        //generate polygon of x4 atoms, and then m and x3 atoms
        Cart3D vx, vy;
        vx=uvPar*parLen;
        vy=Cart3D(uvPer*potentialModel[3],-potentialModel[3]);
        //X4
        int index=0;
        Cart3D crd=atoms[atomIdL].coordinate;
        for(int i=0; i<(nNewX4-1)/2; ++i){
            crd+=vy;
            atoms[x4AtomIds[index]].coordinate=crd;
            ++index;
        }
        crd+=vx*0.5+vy;
        atoms[x4AtomIds[index]].coordinate=crd;
        ++index;
        crd+=vx*0.5-vy;
        for(int i=0; i<nNewX4/2; ++i){
            atoms[x4AtomIds[index]].coordinate=crd;
            ++index;
            crd-=vy;
        }
        //M
        index=0;
        vx=uvPar*(parLen+potentialModel[1]);
        crd=atoms[atomIdL].coordinate-uvPar*potentialModel[1]*0.5+vy*0.5;
        for(int i=0; i<nNewM/2-1; ++i){
            atoms[mAtomIds[index]].coordinate=crd;
            ++index;
            crd+=vy;
        }
        crd+=vy;
        for(int i=0; i<2; ++i){
            crd+=vx/3.0;
            atoms[mAtomIds[index]].coordinate=crd;
            ++index;
        }
        crd+=vx/3.0-vy*2.0;
        for(int i=0; i<nNewM/2-1; ++i){
            atoms[mAtomIds[index]].coordinate=crd;
            ++index;
            crd-=vy;
        }
        //X3
        index=0;
        vx=uvPar*(parLen+potentialModel[1]*3.0);
        crd=atoms[atomIdL].coordinate-uvPar*potentialModel[1]*1.5+vy*0.5;
        for(int i=0; i<nNewX3/2-1; ++i){
            atoms[x3AtomIds[index]].coordinate=crd;
            ++index;
            crd+=vy;
        }
        crd+=vy*2.0;
        for(int i=0; i<2; ++i){
            crd+=vx/3.0;
            atoms[x3AtomIds[index]].coordinate=crd;
            ++index;
        }
        crd+=vx/3.0-vy*3.0;
        for(int i=0; i<nNewX3/2-1; ++i){
            atoms[x3AtomIds[index]].coordinate=crd;
            ++index;
            crd-=vy;
        }
    }
    else if (nNewTriangles%2==1){
        //generate polygon of x4 atoms, and then m and x3 atoms
        int index=0;
        Cart3D vx, vy;
        vx=uvPar*parLen;
        vy=Cart3D(uvPer*potentialModel[3],-potentialModel[3]);
        //X4
        Cart3D crd=atoms[atomIdL].coordinate;
        for(int i=0; i<nNewX4/2; ++i){
            crd+=vy;
            atoms[x4AtomIds[index]].coordinate=crd;
            ++index;
        }
        crd+=vx;
        for(int i=0; i<nNewX4/2; ++i){
            atoms[x4AtomIds[index]].coordinate=crd;
            ++index;
            crd-=vy;
        }
        //M
        index=0;
        vx=uvPar*(parLen+potentialModel[1]);
        crd=atoms[atomIdL].coordinate-uvPar*potentialModel[1]*0.5+vy*0.5;
        for(int i=0; i<(nNewM-1)/2; ++i){
            atoms[mAtomIds[index]].coordinate=crd;
            ++index;
            crd+=vy;
        }
        crd+=vx*0.5;
        atoms[mAtomIds[index]].coordinate=crd;
        ++index;
        crd+=vx*0.5-vy;
        for(int i=0; i<(nNewM-1)/2; ++i){
            atoms[mAtomIds[index]].coordinate=crd;
            ++index;
            crd-=vy;
        }
        //X3
        index=0;
        vx=uvPar*(parLen+potentialModel[1]*3.0);
        crd=atoms[atomIdL].coordinate-uvPar*potentialModel[1]*1.5+vy*0.5;
        for(int i=0; i<(nNewX3-1)/2; ++i){
            atoms[x3AtomIds[index]].coordinate=crd;
            ++index;
            crd+=vy;
        }
        crd+=vx*0.5+vy;
        atoms[x3AtomIds[index]].coordinate=crd;
        ++index;
        crd+=vx*0.5-vy*2.0;
        for(int i=0; i<(nNewM-1)/2; ++i){
            atoms[x3AtomIds[index]].coordinate=crd;
            ++index;
            crd-=vy;
        }
    }
    return;
}
void NetworkCart3DS::popRing(int ringSize, vector<int> &unitPath){
     //remove last built ring

    //calculate number of species
    int nNewTriangles=ringSize-unitPath.size();
    int nNewM=nNewTriangles;
    int nNewX=3*nNewM-1-nNewM;

    //remove atoms, units and ring
    for(int i=0; i<nNewM; ++i) delAtom();
    for(int i=0; i<nNewX; ++i) delAtom();
    for(int i=0; i<nNewTriangles; ++i) delUnit();
    delRing();

    //decrease coordination of dangling atoms in unit path
    int atomIdL = boundaryStatus[find(boundaryUnits.begin(), boundaryUnits.end(), unitPath[0]) - boundaryUnits.begin()];
    int atomIdR = boundaryStatus[find(boundaryUnits.begin(), boundaryUnits.end(), unitPath.rbegin()[0]) - boundaryUnits.begin()];
    --atoms[atomIdL].coordination;
    --atoms[atomIdR].coordination;

    //remove triangle-triangle connections
    delUnitUnitCnx(unitPath[0],nUnits);
    delUnitUnitCnx(unitPath.rbegin()[0],nUnits+nNewTriangles-1);

    //remove triangle-ring connections
    for(int i=0; i<unitPath.size(); ++i) delUnitRingCnx(unitPath[i],nRings);

    //remove ring-ring connections
    vector<int> nbRings;
    nbRings.clear();
    for(int i=0; i<unitPath.size(); ++i){
        for(int j=0; j<units[unitPath[i]].rings.n;++j){
            nbRings.push_back(units[unitPath[i]].rings.ids[j]);
        }
    }
    sort(nbRings.begin(), nbRings.end());
    int prevRing=-1, currRing;
    for(int i=0; i<nbRings.size(); ++i){
        currRing=nbRings[i];
        if(currRing!=prevRing && currRing!=nRings){
            delRingRingCnx(currRing,nRings);
            prevRing=currRing;
        }
    }
    return;
}

void NetworkCart3DS::checkOverlap(){
     //check for overlap of any triangles by projecting onto x-y plane

    unitOverlap=false;

    //make list of all lines that make up triangles
    int a, b, c;
    vector<int> lines;
    lines.clear();
    for(int i=0; i<nUnits; ++i){
        a=units[i].atomsX.ids[0];
        b=units[i].atomsX.ids[1];
        c=units[i].atomsX.ids[2];
        lines.push_back(a);
        lines.push_back(b);
        lines.push_back(a);
        lines.push_back(c);
        lines.push_back(b);
        lines.push_back(c);
    }

    //check overlap of all pairs of lines
    int n=lines.size()/2;
    double x0,x1,x2,x3,y0,y1,y2,y3;
    for(int i=0; i<n-1; ++i){
        x0=atoms[lines[2*i]].coordinate.x;
        y0=atoms[lines[2*i]].coordinate.y;
        x1=atoms[lines[2*i+1]].coordinate.x;
        y1=atoms[lines[2*i+1]].coordinate.y;
        for(int j=i+1; j<n; ++j){
            x2=atoms[lines[2*j]].coordinate.x;
            y2=atoms[lines[2*j]].coordinate.y;
            x3=atoms[lines[2*j+1]].coordinate.x;
            y3=atoms[lines[2*j+1]].coordinate.y;
            if(properIntersectionLines(x0,y0,x1,y1,x2,y2,x3,y3)){
                unitOverlap=true;
                break;
            }
        }
    }
    return;
}

void NetworkCart3DS::calculatePercolation(string shape) {
    //calculate cluster sizes and if percolation occurs
    return;
}