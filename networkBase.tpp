#include "networkBase.h"


//##### NETWORK BASE #####

template <typename CrdT>
Network<CrdT>::Network() {
    //default constructor
    nAtoms=0;
    nUnits=0;
    nRings=0;
    energy=numeric_limits<double>::infinity();
    optIterations=-1;
    atoms.clear();
    units.clear();
    rings.clear();
    boundaryUnits.clear();
}

template <typename CrdT>
int Network<CrdT>::getNRings() {
    //return number of rings in network
    return nRings;
}

template <typename CrdT>
double Network<CrdT>::getEnergy() {
    //return energy of network
    return energy;
}

template <typename CrdT>
int Network<CrdT>::getIterations() {
    //return minimisation iterations
    return optIterations;
}

template <typename CrdT>
void Network<CrdT>::addAtom(Atom<CrdT> atom) {
    //add atom to network and update map
    atoms.push_back(atom);
    ++nAtoms;
}

template <typename CrdT>
void Network<CrdT>::addUnit(Unit unit) {
    //add unit to network and update map
    units.push_back(unit);
    ++nUnits;
}

template <typename CrdT>
void Network<CrdT>::addRing(Ring ring) {
    //add ring to network and update map
    rings.push_back(ring);
    ++nRings;
}

template <typename CrdT>
void Network<CrdT>::delAtom() {
    //delete atom from network
    atoms.pop_back();
    --nAtoms;
}

template <typename CrdT>
void Network<CrdT>::delUnit() {
    //delete unit from network
    units.pop_back();
    --nUnits;
}

template <typename CrdT>
void Network<CrdT>::delRing() {
    //delete ring from network
    rings.pop_back();
    --nRings;
}

template <typename CrdT>
int Network<CrdT>::addUnitAtomXCnx(int uId, int aId) {
    //add atom x id to unit
    int status=units[uId].atomsX.add(aId);
    return status;
}

template <typename CrdT>
int Network<CrdT>::addUnitRingCnx(int uId, int rId) {
    //add mutual unit-ring connection
    int uStatus=units[uId].rings.add(rId);
    int rStatus=rings[rId].units.add(uId);
    return uStatus+rStatus;
}

template <typename CrdT>
int Network<CrdT>::addUnitUnitCnx(int uId1, int uId2) {
    //add mutual unit-unit connection
    int status1=units[uId1].units.add(uId2);
    int status2=units[uId2].units.add(uId1);
    return status1+status2;
}

template <typename CrdT>
int Network<CrdT>::addRingRingCnx(int rId1, int rId2) {
    //add mutual unit-unit connection
    int status1=rings[rId1].rings.add(rId2);
    int status2=rings[rId2].rings.add(rId1);
    return status1+status2;
}

template <typename CrdT>
void Network<CrdT>::delUnitUnitCnx(int uId1, int uId2) {
    //del unit-unit connection
    units[uId1].units.del(uId2);
}

template <typename CrdT>
void Network<CrdT>::delUnitRingCnx(int uId, int rId) {
    //delete unit-ring connection
    units[uId].rings.del(rId);
}

template <typename CrdT>
void Network<CrdT>::delRingRingCnx(int rId1, int rId2) {
    //delete ring-ring connection
    rings[rId1].rings.del(rId2);
}

template <typename CrdT>
void Network<CrdT>::changeUnitAtomXCnx(int uId, int aId1, int aId2) {
    //change unit-atom x connection
    units[uId].atomsX.change(aId1,aId2);
}


template <typename CrdT>
bool Network<CrdT>::checkActiveUnit(int &uId, int sumCheck) {
    //checks if a unit is active by summing coordination of associated x atoms
    double cndSum=0;
    for(int i=0; i<units[uId].atomsX.n; ++i){
        cndSum+=atoms[units[uId].atomsX.ids[i]].coordination;
    }
    if(cndSum<sumCheck) return true;
    else return false;
}

template <typename CrdT>
bool Network<CrdT>::checkEdgeUnit(int &uId, int ringCheck) {
    //checks if a unit is on edge by checking number of associated rings
    if(units[uId].rings.n<ringCheck) return true;
    else return false;
}

template <typename CrdT>
bool Network<CrdT>::trialRing(int ringSize, vector<int> &unitPath, vector<double> &potentialModel) {
    //build a ring of a given size to a starting path
    //minimise and calculate energy
    //remove ring

    //check build type - path closure or need to add units
    bool ring0; //flag for path closure, auxiliary index
    if(ringSize==unitPath.size()) ring0=true;
    else ring0=false;

    //build trial ring
    if(ring0) buildRing0(unitPath);
    else buildRing(ringSize, unitPath, potentialModel);

    //geometry optimise
    geometryOptimiseLocal(potentialModel);
//    geometryOptimiseGlobal(potentialModel);

    //check for geometry anomalies
    bool geometryCheck=checkLocalGrowth(nRings-1);

    //pop trial ring
    if(ring0) popRing0(unitPath);
    else popRing(ringSize, unitPath);

    return geometryCheck;
}

template <typename CrdT>
void Network<CrdT>::acceptRing(int ringSize, vector<int> &unitPath, vector<double> &potentialModel) {
    //build ring of given size to a starting path, minimise and calculate boundary

    if(ringSize==unitPath.size()) buildRing0(unitPath);
    else buildRing(ringSize, unitPath, potentialModel);
    geometryOptimiseLocal(potentialModel);
//    geometryOptimiseGlobal(potentialModel);
    calculateBoundary();
}

template <typename CrdT>
void Network<CrdT>::clean() {
    //remove any uncoordinated atoms and reassign ids

    //loop over atoms and find zero-coordinate species
    map<int,int> updatedAtomIds; //reassigned ids
    vector<int> removeAtoms; //atoms with zero coordination - artefacts of build process from buildRing0
    removeAtoms.clear();
    int aId=0;
    for(int i=0; i<nAtoms; ++i){
        if(atoms[i].coordination==0){
            removeAtoms.push_back(i);
            updatedAtomIds[i]=-1;
        }
        else{
            updatedAtomIds[i]=aId;
            ++aId;
        }
    }

    //remove atoms from list in reverse order
    for(int i=0; i<removeAtoms.size(); ++i) atoms.erase(atoms.begin()+removeAtoms.rbegin()[i]);
    nAtoms-=removeAtoms.size();

    //update atom ids in unit connections
    int aId0, aId1;
    for(int i=0; i<nUnits; ++i){
        for(int j=0; j<units[i].atomsX.n; ++j){//update x atom ids
            aId0=units[i].atomsX.ids[j];
            aId1=updatedAtomIds[aId0];
            if(aId1==-1) cout<<"ERROR IN NETWORK CLEANING"<<endl;
            units[i].atomsX.ids[j]=aId1;
        }
        //update m atom id
        aId0=units[i].atomM;
        aId1=updatedAtomIds[aId0];
        if(aId1==-1) cout<<"ERROR IN NETWORK CLEANING"<<endl;
        units[i].atomM=aId1;
    }
}

template <typename CrdT>
void Network<CrdT>::calculateBoundary() {
    //find units on boundary of network
    //follow units on edge, flag when traversed
    //an edge triangle will have rings.n<3
    //an active triangle will have an atom with coordination<4
    //if active must be on edge

    boundaryUnits.clear();
    boundaryStatus.clear();

    //starting position as an active triangle on edge and arbitrary edge neighbour
    for(int i=0; i<nUnits; ++i){
        if(checkActiveUnit(i)){
            boundaryUnits.push_back(i);
            boundaryStatus.push_back(true);
            units[i].flag=true;
            //must have two neighbours on edge so pick first arbitrarily
            boundaryUnits.push_back(units[i].units.ids[0]);
            boundaryStatus.push_back(checkActiveUnit(boundaryUnits[1]));
            units[boundaryUnits.rbegin()[0]].flag=true;
            break;
        }
    }

    //trace perimeter back to start
    int id0=boundaryUnits.rbegin()[0], id1, id1a, id1b, id1c, id1d, id1e;
    bool active0=boundaryStatus.rbegin()[0], active1;
    bool complete=false;
    for(;;){
        if(id0==16){
            int a=0;
        }
        if(active0){//only one possiblity - pick non-flagged path
            for(int i=0; i<2; ++i){
                id1=units[id0].units.ids[i];
                if(!units[id1].flag) break;
                else if(id1==boundaryUnits[0] && boundaryUnits.size()>2) complete=true;
            }
        }
        else{
            //get ids of three connected units
            id1c=units[id0].units.ids[0];
            id1d=units[id0].units.ids[1];
            id1e=units[id0].units.ids[2];
            if(id1c==boundaryUnits[0] && boundaryUnits.size()>2) complete=true;
            if(id1d==boundaryUnits[0] && boundaryUnits.size()>2) complete=true;
            if(id1e==boundaryUnits[0] && boundaryUnits.size()>2) complete=true;
            //keep ids of those which are not flagged
            id1a=-1;
            id1b=-1;
            if(!units[id1c].flag){
                id1a=id1c;
                if(!units[id1d].flag) id1b=id1d;
                else if(!units[id1e].flag) id1b=id1e;
            }
            else if(!units[id1d].flag){
                id1a=id1d;
                if(!units[id1e].flag) id1b=id1e;
            }
            else if(!units[id1e].flag) id1a=id1e;
            //review possiblities
            if(id1b==-1){//only one non-flag option
                id1=id1a;
            }
            else{//two non-flag options, find number of edge options
                bool edgeA=checkEdgeUnit(id1a);
                bool edgeB=checkEdgeUnit(id1b);
                if(edgeA && !edgeB){//only a on edge
                    id1=id1a;
                }
                else if(!edgeA && edgeB){//only b on edge
                    id1=id1b;
                }
                else{//two edge options
                    //check number of associated rings - pick unit that has just one
                    if(units[id1a].rings.n==1 && units[id1b].rings.n>1) id1=id1a;
                    else if(units[id1a].rings.n>1 && units[id1b].rings.n==1) id1=id1b;
                    else{
                        //both have multiple rings - hardest case but also most rare
                        //compare which rings each unit belongs to, pick unit with just one match
                        int ring0i=units[id0].rings.ids[0];
                        int ring0ii=units[id0].rings.ids[1];
                        int ring1ai=units[id1a].rings.ids[0];
                        int ring1aii=units[id1a].rings.ids[1];
                        int ring1bi=units[id1b].rings.ids[0];
                        int ring1bii=units[id1b].rings.ids[1];
                        if(ring1ai!=ring0i && ring1ai!=ring0ii) id1=id1a;
                        else if(ring1aii!=ring0i && ring1aii!=ring0ii) id1=id1a;
                        else if(ring1bi!=ring0i && ring1bi!=ring0ii) id1=id1b;
                        else if(ring1bii!=ring0i && ring1bii!=ring0ii) id1=id1b;
                        else{
                            cout<<"Boundary failed"<<endl;
//                            consoleVector(boundaryUnits);
                            Logfile dump;
                            write("dump",false,dump);
                            exit(9);
                        }
                    }
                }
            }
        }
        if (complete) break;
        else {
            active1=checkActiveUnit(id1);
            boundaryUnits.push_back(id1);
            boundaryStatus.push_back(active1);
            units[id1].flag=true;
            id0 = id1;
            active0=active1;
        }
    }

    //convert boundary status values to atom which is active
    for(int i=0; i<boundaryStatus.size(); ++i){
        if(boundaryStatus[i]==0) boundaryStatus[i]=-1; //as 0 can be an id
        else{
            id0=boundaryUnits[i];
            for(int j=0; j<3; ++j){//loop over x atoms and find undercoordinated atom
                id1=units[id0].atomsX.ids[j];
                if(atoms[id1].coordination<4){
                    boundaryStatus[i]=id1;
                    break;
                }
            }
        }

    }

    //remove flag
    for(int i=0; i<boundaryUnits.size(); ++i) units[boundaryUnits[i]].flag=false;

//    consoleVector(boundaryUnits);
//    consoleVector(boundaryStatus);
}

template <typename CrdT>
vector<int> Network<CrdT>::getBoundarySection(int startId, bool direction) {
    //find section of unit boundary in given direction

    //find position of id in boundary
    int startPos = find(boundaryUnits.begin(), boundaryUnits.end(), startId) - boundaryUnits.begin();

    vector<int> section;
    section.clear();
    section.push_back(startId);
    int n=boundaryUnits.size();

    //search in one of two directions, loop over perimeter until find next active unit
    if(direction) {
        int j;
        for (int i = 1; i <= n; ++i) {
            j=(startPos + i) % n;
            section.push_back(boundaryUnits[j]);
            if(boundaryStatus[j]>=0) break;
        }
    }
    else{
        int j;
        for (int i = 1; i <= n; ++i) {
            j = (startPos + n - i) % n;
            section.push_back(boundaryUnits[j]);
            if(boundaryStatus[j]>=0) break;
        }
    }
    return section;
}

template <typename CrdT>
void Network<CrdT>::findLocalRegion(int &rId, int nFlexShells) {
    //find units around and included in a given ring, within a given number of connections, make map of corresponding local atoms
    //have flexible shells, then fixed shell

    localAtomMap.clear();
    globalAtomMap.clear();
    flexLocalUnits.clear();
    fixedLocalUnits.clear();
    fixedLocalAtoms.clear();

    //get units in ring
    vector<int> shell0(rings[rId].units.n), shell1;
    for(int i=0; i<rings[rId].units.n; ++i) shell0[i]=rings[rId].units.ids[i];
    for(int i=0; i<shell0.size(); ++i) flexLocalUnits.push_back(shell0[i]);

    //loop over flexible shells and get units, then get fixed shell
    for(int i=0; i<nFlexShells+1; ++i){
        shell1.clear();
        //find adjecent units to shell0
        for(int j=0; j<shell0.size(); ++j){
            for(int k=0; k<units[shell0[j]].units.n; ++k){
                shell1.push_back(units[shell0[j]].units.ids[k]);
            }
        }
        if(shell1.size()>0){
            //get unique units not in shell0
            sort(shell1.begin(), shell1.end());
            shell1.erase(unique(shell1.begin(), shell1.end()),shell1.end());
            for (int i=0; i <flexLocalUnits.size(); ++i) shell1.erase(remove(shell1.begin(), shell1.end(), flexLocalUnits[i]), shell1.end());
            //add to vector
            if(i!=nFlexShells){
                for(int i=0; i<shell1.size(); ++i) flexLocalUnits.push_back(shell1[i]);
            }
            else{
                for(int i=0; i<shell1.size(); ++i) fixedLocalUnits.push_back(shell1[i]);
            }
            shell0=shell1;
        }
        else break;
    }

    //make map of atoms to include in local region
    int m, x;
    nLocalAtoms=0;
    for(int i=0; i<flexLocalUnits.size(); ++i){
        m=units[flexLocalUnits[i]].atomM;
        if(localAtomMap.count(m)==0){
            localAtomMap[m]=nLocalAtoms;
            globalAtomMap[nLocalAtoms]=m;
            ++nLocalAtoms;
        }
        for(int j=0; j<units[flexLocalUnits[i]].atomsX.n; ++j){
            x=units[flexLocalUnits[i]].atomsX.ids[j];
            if(localAtomMap.count(x)==0){
                localAtomMap[x]=nLocalAtoms;
                globalAtomMap[nLocalAtoms]=x;
                ++nLocalAtoms;
            }
        }
    }
    for(int i=0; i<fixedLocalUnits.size(); ++i){
        m=units[fixedLocalUnits[i]].atomM;
        if(localAtomMap.count(m)==0){
            localAtomMap[m]=nLocalAtoms;
            globalAtomMap[nLocalAtoms]=m;
            fixedLocalAtoms.push_back(nLocalAtoms);
            ++nLocalAtoms;
        }
        for(int j=0; j<units[fixedLocalUnits[i]].atomsX.n; ++j){
            x=units[fixedLocalUnits[i]].atomsX.ids[j];
            if(localAtomMap.count(x)==0){
                localAtomMap[x]=nLocalAtoms;
                globalAtomMap[nLocalAtoms]=x;
                fixedLocalAtoms.push_back(nLocalAtoms);
                ++nLocalAtoms;
            }
        }
    }

//    cout<<"***"<<endl;
//    for(int i=0; i<fixedLocalAtoms.size(); ++i) {
//        cout<<globalAtomMap[fixedLocalAtoms[i]]<<endl;
//    }
//    cout<<"***"<<endl;
//
//    cout<<"-----"<<endl;
//    for(int i=0; i<nLocalAtoms;++i){
//        cout<<globalAtomMap[i]<<endl;
//    }
//    cout<<"-----"<<endl;


}

template <typename CrdT>
void Network<CrdT>::calculateRingStatistics() {
    //calculate ring statistics, ring statistics around each ring, and aboav-weaire analysis

    //calculate distribution of ring sizes and store unique ring sizes
    vector<int> ringSizes;
    for(int i=0; i<nRings; ++i) ringSizes.push_back(rings[i].units.n);
    DiscreteDistribution ringStats(ringSizes);
    ringStatistics=ringStats;

    //calculate distribution of ring sizes excluding edge rings
    ringSizes.clear();
    for(int i=0; i<nRings; ++i){
        if(rings[i].rings.full) ringSizes.push_back(rings[i].units.n);
    }
    ringStats=DiscreteDistribution(ringSizes);
    bulkRingStatistics=ringStats;

    //calculate distributions for each ring size (excluding edge rings)
    int ringRef;
    indRingStatistics.clear();
    for(int i=0; i<ringStatistics.n; ++i){//loop over ring sizes
        int s=ringStatistics.x[i];
        ringSizes.clear();
        for(int j=0; j<nRings; ++j){//get ring sizes around ring of given size
            if(rings[j].rings.full){//only include rings not on edge
                if(rings[j].units.n==s){
                    for(int k=0; k<rings[j].rings.n; ++k){
                        ringSizes.push_back(rings[rings[j].rings.ids[k]].units.n);
                    }
                }
            }
        }
        if(ringSizes.size()>0){
            DiscreteDistribution ringStats(ringSizes);
            indRingStatistics[s]=ringStats;
        }
    }

    //calculate aboav-weaire fit
    vector<double> x;
    vector<double> y;
    for(int i=0; i<ringStatistics.n; ++i){
        int s=ringStatistics.x[i];
        if(indRingStatistics.count(s)>0){
            x.push_back(ringStatistics.mean*(s-ringStatistics.mean));
            y.push_back(s*indRingStatistics[s].mean);
        }
    }
    aboavWeaireParameters=leastSquaresLinearRegression(x,y);
    aboavWeaireParameters[0]=1.0-aboavWeaireParameters[0]; //alpha
    aboavWeaireParameters[1]-=ringStatistics.mean*ringStatistics.mean; //mu

    //initialise ring colours
    ringColours.resize(nRings,col_vector<int>(4));
    for(int i=0; i<nRings; ++i) ringColours[i][0]=rings[i].units.n;
}

template <typename CrdT>
void Network<CrdT>::calculateRingAreas() {
    //calculate dimensionless areas of rings separated by size

    //loop over rings, calculate dimensionless area and store in vector according to ring size
    map<int, vector<double> > ringSizeAreas;
    double area;
    double mm_sq=bondLenDistMM.mean*bondLenDistMM.mean; //to make dimensionless
    int ringSize,m0, m1;
    for(int i=0; i<nRings; ++i){
        area=0.0;
        ringSize=rings[i].units.n;
        for(int j=0; j<ringSize; ++j){
            m0=units[rings[i].units.ids[j]].atomM;
            m1=units[rings[i].units.ids[(j+1)%ringSize]].atomM;
            area+=atoms[m0].coordinate.x*atoms[m1].coordinate.y;
            area-=atoms[m1].coordinate.x*atoms[m0].coordinate.y;
        }
        area=fabs(0.5*area/mm_sq);
        ringSizeAreas[ringSize].push_back(area);
    }

    //make distributions according to ring size
    vector<int> ringSizes=ringStatistics.getValues();
    for(int i=0; i<ringSizes.size(); ++i){
        ContinuousDistribution areaDistribution(ringSizeAreas.at(ringSizes[i]));
        ringAreas[ringSizes[i]]=areaDistribution;
    }

}

template <typename CrdT>
void Network<CrdT>::calculateBondDistributions(bool fullDist) {
    //calculate bond length/angle distributions

    //flag whether to write full distributions
    writeFullDistributions=fullDist;

    //M-X length
    vector<double> bondLengths;
    bondLengths.clear();
    int m,x;
    CrdT crdM, crdX, crdMX;
    for(int i=0; i<nUnits; ++i){
        m=units[i].atomM;
        crdM=atoms[m].coordinate;
        for(int j=0; j<units[i].atomsX.n; ++j){
            x=units[i].atomsX.ids[j];
            crdX=atoms[x].coordinate;
            crdMX=crdM-crdX;
            bondLengths.push_back(crdMX.norm());
        }
    }
    ContinuousDistribution bondLenMX(bondLengths);
    bondLenDistMX=bondLenMX;

    //X-X length
    bondLengths.clear();
    int x0, x1;
    CrdT crdX0, crdX1, crdXX;
    for(int i=0; i<nUnits; ++i){
        for(int j=0; j<units[i].atomsX.n-1; ++j){
            x0=units[i].atomsX.ids[j];
            crdX0=atoms[x0].coordinate;
            for(int k=j+1; k<units[i].atomsX.n; ++k){
                x1=units[i].atomsX.ids[k];
                crdX1=atoms[x1].coordinate;
                crdXX=crdX1-crdX0;
                bondLengths.push_back(crdXX.norm());
            }
        }
    }
    ContinuousDistribution bondLenXX(bondLengths);
    bondLenDistXX=bondLenXX;

    //M-M length
    bondLengths.clear();
    int m0, m1;
    CrdT crdM0, crdM1, crdMM;
    for(int i=0; i<nUnits; ++i){
        m0=units[i].atomM;
        crdM0=atoms[m0].coordinate;
        for(int j=0; j<units[i].units.n; ++j){
            m1=units[units[i].units.ids[j]].atomM;
            if(m0<m1){//prevent double  counting
                crdM1=atoms[m1].coordinate;
                crdMM=crdM1-crdM0;
                bondLengths.push_back(crdMM.norm());
            }
        }
    }
    ContinuousDistribution bondLenMM(bondLengths);
    bondLenDistMM=bondLenMM;

    //M-X-M angle
    vector<double> bondAngles;
    bondAngles.clear();
    int u0,u1;
    int x00, x01, x02, x10, x11, x12;
    CrdT crdMX0, crdMX1;
    double theta0, theta1;
    for(int i=0; i<nUnits; ++i){
        u0=i;
        m0=units[u0].atomM;
        x00=units[u0].atomsX.ids[0];
        x01=units[u0].atomsX.ids[1];
        x02=units[u0].atomsX.ids[2];
        crdM0=atoms[m0].coordinate;
        for(int j=0; j<units[i].units.n; ++j){
            u1=units[i].units.ids[j];
            if(u0<u1){//prevent double counting
                m1=units[u1].atomM;
                x10=units[u1].atomsX.ids[0];
                x11=units[u1].atomsX.ids[1];
                x12=units[u1].atomsX.ids[2];
                crdM1=atoms[m1].coordinate;
                //find bridging x atom
                if(x00==x10) x=x00;
                else if(x01==x10) x=x01;
                else if(x02==x10) x=x02;
                else if(x00==x11) x=x00;
                else if(x01==x11) x=x01;
                else if(x02==x11) x=x02;
                else if(x00==x12) x=x00;
                else if(x01==x12) x=x01;
                else if(x02==x12) x=x02;
                else cout<<"ERROR IN BOND ANGLE CALCULATION"<<endl;
                crdX=atoms[x].coordinate;
                crdMX0=crdM0-crdX;
                crdMX1=crdM1-crdX;
                crdMX0.normalise();
                crdMX1.normalise();
                theta0=acos(crdMX0*crdMX1);
                theta1=2.0*M_PI-theta0;
                bondAngles.push_back(theta0);
                bondAngles.push_back(theta1);
            }
        }
    }
    ContinuousDistribution bondAngMXM(bondAngles);
    bondAngDistMXM=bondAngMXM;

    //M-M-M angle
    bondAngles.clear();
    double theta;
    CrdT crdM2, crdMM0, crdMM1, crdMM2;
    for(int i=0; i<nUnits; ++i){
        crdM=atoms[units[i].atomM].coordinate;
        if(units[i].units.n==2){
            crdM0=atoms[units[units[i].units.ids[0]].atomM].coordinate;
            crdM1=atoms[units[units[i].units.ids[1]].atomM].coordinate;
            crdMM0=crdM0-crdM;
            crdMM1=crdM1-crdM;
            crdMM0.normalise();
            crdMM1.normalise();
            theta=acos(crdMM0*crdMM1);
            bondAngles.push_back(theta);
        }
        else if(units[i].units.n==3){
            crdM0=atoms[units[units[i].units.ids[0]].atomM].coordinate;
            crdM1=atoms[units[units[i].units.ids[1]].atomM].coordinate;
            crdM2=atoms[units[units[i].units.ids[2]].atomM].coordinate;
            crdMM0=crdM0-crdM;
            crdMM1=crdM1-crdM;
            crdMM2=crdM2-crdM;
            crdMM0.normalise();
            crdMM1.normalise();
            crdMM2.normalise();
            theta=acos(crdMM0*crdMM1);
            bondAngles.push_back(theta);
            theta=acos(crdMM1*crdMM2);
            bondAngles.push_back(theta);
            theta=acos(crdMM0*crdMM2);
            bondAngles.push_back(theta);
        }
        else cout<<"ERROR IN BOND ANGLE CALCULATION"<<endl;
    }
    ContinuousDistribution bondAngMMM(bondAngles);
    bondAngDistMMM=bondAngMMM;

}

template <typename CrdT>
void Network<CrdT>::write(string prefix, bool special, Logfile &logfile) {
    //write network and analysis to files
    writeNetwork(prefix,logfile);
    if(special) writeNetworkSpecial(prefix,logfile);
    writeAnalysis(prefix,logfile);
}

template <typename CrdT>
void Network<CrdT>::writeAnalysis(string prefix, Logfile &logfile) {
    //write analysis to file

    logfile.log("Writing network analysis","","",0,false);
    //set up analysis file
    string analysisFilename = prefix + "_analysis.out";
    ofstream analysisFile(analysisFilename,ios::in|ios::trunc);

    //overlap
    writeFileValue(analysisFile,"Valid geometry",true);
    writeFileValue(analysisFile,!unitOverlap,true);

    //ring statistics
    writeFileValue(analysisFile,"p_n, <n>, s for total, bulk, individual",true);
    analysisFile << fixed << showpoint << setprecision(1);
    //total ring statistics
    vector<int> ringSizes=ringStatistics.getValues();
    writeFileVector(analysisFile,ringSizes);
    analysisFile << fixed << showpoint << setprecision(6);
    vector<double> data=ringStatistics.getProbabilities();
    data.push_back(ringStatistics.mean);
    data.push_back(double(ringStatistics.sampleSize));
    writeFileVector(analysisFile,data);
    //bulk rings statistics
    data.clear();
    for(int i=0; i<ringSizes.size(); ++i) data.push_back(bulkRingStatistics.getProbability(ringSizes[i]));
    data.push_back(bulkRingStatistics.mean);
    data.push_back(double(bulkRingStatistics.sampleSize));
    writeFileVector(analysisFile,data);
    //individual ring statistics
    for(int i=0; i<ringSizes.size(); ++i){
        data.clear();
        int s=ringSizes[i];
        if(indRingStatistics.count(s)==0){
            for(int j=0; j<ringSizes.size()+2; ++j) data.push_back(0.0);
        }
        else{
            for(int j=0; j<ringSizes.size(); ++j){
                data.push_back(indRingStatistics.at(s).getProbability(ringSizes[j]));
            }
            data.push_back(indRingStatistics.at(s).mean);
            data.push_back(double(indRingStatistics.at(s).sampleSize));
        }
        writeFileVector(analysisFile,data);
    }
    logfile.log("Ring statistics written to: ",analysisFilename,"", 1, false);

    //ring correlations
    writeFileValue(analysisFile,"Aboav-Weaire alpha, mu and rsq",true);
    data=aboavWeaireParameters;
    writeFileVector(analysisFile,data);
    logfile.log("Aboav-Weaire parameters written to: ",analysisFilename,"", 1, false);

    //bond length/angle distributions
    writeFileValue(analysisFile,"Bond Length/Angle Distributions: MX, XX, MM, MXM, MMM",true);
    writeFileValue(analysisFile,bondLenDistMX.mean,bondLenDistMX.sdev);
    writeFileValue(analysisFile,bondLenDistXX.mean,bondLenDistXX.sdev);
    writeFileValue(analysisFile,bondLenDistMM.mean,bondLenDistMM.sdev);
    writeFileValue(analysisFile,bondAngDistMXM.mean,bondAngDistMXM.sdev);
    writeFileValue(analysisFile,bondAngDistMMM.mean,bondAngDistMMM.sdev);
    logfile.log("Bond length and angle distribution summaries written to: ",analysisFilename,"",1,false);
    if(writeFullDistributions){
        string mxAnalysisFilename = prefix + "_analysis_mx.out";
        string xxAnalysisFilename = prefix + "_analysis_xx.out";
        string mmAnalysisFilename = prefix + "_analysis_mm.out";
        string mxmAnalysisFilename = prefix + "_analysis_mxm.out";
        string mmmAnalysisFilename = prefix + "_analysis_mmm.out";
        ofstream mxAnalysisFile(mxAnalysisFilename,ios::in|ios::trunc);
        ofstream xxAnalysisFile(xxAnalysisFilename,ios::in|ios::trunc);
        ofstream mmAnalysisFile(mmAnalysisFilename,ios::in|ios::trunc);
        ofstream mxmAnalysisFile(mxmAnalysisFilename,ios::in|ios::trunc);
        ofstream mmmAnalysisFile(mmmAnalysisFilename,ios::in|ios::trunc);
        vector<double> rawDistribution=bondLenDistMX.getValues();
        writeFileVectorTranspose(mxAnalysisFile,rawDistribution);
        rawDistribution=bondLenDistXX.getValues();
        writeFileVectorTranspose(xxAnalysisFile,rawDistribution);
        rawDistribution=bondLenDistMM.getValues();
        writeFileVectorTranspose(mmAnalysisFile,rawDistribution);
        rawDistribution=bondAngDistMXM.getValues();
        writeFileVectorTranspose(mxmAnalysisFile,rawDistribution);
        rawDistribution=bondAngDistMMM.getValues();
        writeFileVectorTranspose(mmmAnalysisFile,rawDistribution);
        logfile.log("MX bond length distribution written to: ",mxAnalysisFilename,"",1,false);
        logfile.log("XX bond length distribution written to: ",xxAnalysisFilename,"",1,false);
        logfile.log("MM bond length distribution written to: ",xxAnalysisFilename,"",1,false);
        logfile.log("MXM bond angle distribution written to: ",mxmAnalysisFilename,"",1,false);
        logfile.log("MMM bond angle distribution written to: ",mmmAnalysisFilename,"",1,false);
        mxAnalysisFile.close();
        xxAnalysisFile.close();
        mxmAnalysisFile.close();
    }

    //ring areas
    if(ringAreas.size()>0){//only write if performed analysis
        writeFileValue(analysisFile,"Average Ring Areas",true);
        writeFileVector(analysisFile,ringSizes);
        vector<double> meanRingAreas;
        for(int i=0; i<ringSizes.size(); ++i) meanRingAreas.push_back(ringAreas.at(ringSizes[i]).mean);
        writeFileVector(analysisFile,meanRingAreas);
    }

    //clustering and percolation
    if(clusterDistributions.size()>0){//only write if performed analysis
        writeFileValue(analysisFile,"Raw Cluster Distributions",true);
        for(int i=0; i<ringSizes.size(); ++i){
            writeFileVector(analysisFile,clusterDistributions.at(ringSizes[i]).getValues());
            writeFileVector(analysisFile,clusterDistributions.at(ringSizes[i]).getRawProbabilities());
        }
        writeFileValue(analysisFile,"Percolation",true);
        for(int i=0;i<ringSizes.size();++i) writeFileValue(analysisFile,percolation.at(ringSizes[i]),false);
        writeFileValue(analysisFile,"  ",true);
        logfile.log("Cluster distributions and percolation written to: ",analysisFilename,"", 1, false);
    }

    logfile.log("Writing complete","","",0,true);
}

template <typename CrdT>
void Network<CrdT>::kill(string prefix, Logfile &logfile) {
    //kill network growth early and write out files

    //clean, initialise and get ring colours
    clean();
    ringColours.resize(nRings,col_vector<int>(2));
    for(int i=0; i<nRings; ++i) ringColours[i][0]=rings[i].units.n;

    //write network only - not analysis
    writeNetwork(prefix,logfile);
}
//template <typename CrdT>
//Network<CrdT>::
