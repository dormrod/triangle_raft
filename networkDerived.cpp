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

int NetworkCart2D::getActiveUnit(string shape, double size) {
    //find active unit within given shape, return -1 if cannot find

    int unitId=-1;
    int atomId; //id of atom with dangling bond

    if(shape=="C") {//circle
        double rsq = size * size;
        for (int i = 0; i < boundaryUnits.size(); ++i) {
            atomId = boundaryStatus[i];
            if (atomId >= 0) {
                if (atoms[atomId].coordinate.normSq() < rsq) {
                    unitId = boundaryUnits[i];
                    break;
                }
            }
        }
    }
    else if(shape=="S"){
        Cart2D c;
        for (int i = 0; i < boundaryUnits.size(); ++i) {
            atomId = boundaryStatus[i];
            if (atomId >= 0) {
                c=atoms[atomId].coordinate;
                if (c.x < size && c.x > -size && c.y < size && c.y > - size) {
                    unitId = boundaryUnits[i];
                    break;
                }
            }
        }
    }
    else if(shape=="H"){
        Cart2D c;
        double f;
        for (int i = 0; i < boundaryUnits.size(); ++i) {
            atomId = boundaryStatus[i];
            if (atomId >= 0) {
                c=atoms[atomId].coordinate;
                c/=size;
                f=c.x*c.x+pow((1.25*c.y-sqrt(fabs(c.x))),2);
                if (f<1) {
                    unitId = boundaryUnits[i];
                    break;
                }
            }
        }
    }
    return unitId;
}

void NetworkCart2D::trialRing(int ringSize, vector<int> &unitPath, vector<double> &potentialModel) {
    //build a ring of a given size to a starting path
    //minimise and calculate energy
    //remove ring

    //build trial ring
    buildRing(ringSize, unitPath, potentialModel);

    //geometry optimise
    geometryOptimise(potentialModel);

    //pop trial ring
    popRing(ringSize, unitPath);

}

void NetworkCart2D::popRing(int ringSize, vector<int> &unitPath) {
    //remove last built ring

    //calculate number of species
    int nNewTriangles=ringSize-unitPath.size();
    int nNewM=nNewTriangles;
    int nNewX=3*nNewM-1-nNewM;
    int nNewX3=nNewTriangles;
    int nNewX4=nNewX-nNewX3;

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
}

void NetworkCart2D::buildRing(int ringSize, vector<int> &unitPath, vector<double> &potentialModel) {
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
        Atom<Cart2D> atom(nAtoms,14,3);
        addAtom(atom);
    }
    for(int i=0; i<nNewX3; ++i){//x3
        Atom<Cart2D> atom(nAtoms,8,3);
        addAtom(atom);
    }
    for(int i=0; i<nNewX4; ++i){//x4
        Atom<Cart2D> atom(nAtoms,8,4);
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

    //generate new coordinates
    //generate vectors
    Cart2D uvPar, uvPer; //unit vectors parallel and perpendicular to L->R
    uvPar=atoms[atomIdR].coordinate-atoms[atomIdL].coordinate;
    Cart2D vClockwise(uvPar.y,-uvPar.x);
    Cart2D vAntiClockwise(-uvPar.y,uvPar.x);
    Cart2D dir=atoms[atomIdL].coordinate-atoms[units[unitPath[0]].atomM].coordinate;
    if(vClockwise*dir>0) uvPer=vClockwise;
    else uvPer=vAntiClockwise;
    double perLen, parLen; //lengths of original perpendicular and parallel vectors
    uvPar.normalise(parLen);
    uvPer.normalise(perLen);

    if (nNewTriangles%2==0){
        //generate polygon of x4 atoms, and then m and x3 atoms
        Cart2D vx, vy;
        vx=uvPar*parLen;
        vy=uvPer*potentialModel[3];
        //X4
        int index=0;
        Cart2D crd=atoms[atomIdL].coordinate;
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
        Cart2D vx, vy;
        vx=uvPar*parLen;
        vy=uvPer*potentialModel[3];
        //X4
        Cart2D crd=atoms[atomIdL].coordinate;
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
