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
        if(active0){//only one possiblity - pick non-flagged path
            for(int i=0; i<2; ++i){
                id1=units[id0].units.ids[i];
                active1=checkActiveUnit(id1);
                if(!units[id1].flag) break;
            }
            if(id1==boundaryUnits[0]) complete=true;
        }
        else{
            //get ids of three connected units
            id1c=units[id0].units.ids[0];
            id1d=units[id0].units.ids[1];
            id1e=units[id0].units.ids[2];
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
                            consoleVector(boundaryUnits);
                            Logfile dump;
                            write("dump",dump);
                            exit(9);
                        }
                    }
                }
            }
        }
        if (complete) break;
        else {
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

    consoleVector(boundaryUnits);
    consoleVector(boundaryStatus);

}

