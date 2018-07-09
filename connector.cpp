#include "connector.h"

//##### CONNECTOR #####
Connector::Connector() {
    //default constructor
    n=0;
    max=0;
}

Connector::Connector(int maxN){
    //constructor
    n=0;
    full=false;
    max=maxN;
    ids=new int[max]();
}

Connector::~Connector() {
    //destructor
    if(max>0) delete[] ids;
}

Connector::Connector(const Connector &source) {
    //copy constructor

    //shallow copies
    n=source.n;
    full=source.full;
    max=source.max;

    //deep copies
    if(max>0){
        ids=new int[max]();
        for(int i=0; i<n; ++i) ids[i]=source.ids[i];
    }
}

Connector& Connector::operator=(const Connector &source) {
    //overload assignment operator

    if (this == &source) return *this;

    //shallow copies
    n=source.n;
    full=source.full;
    max=source.max;

    //deep copies
    if(max>0){
        ids=new int[max]();
        for(int i=0; i<n; ++i) ids[i]=source.ids[i];
    }

    return *this;
}

int Connector::add(int cnx) {
    //add a connection if not full, check if full after addition

    if(full) return 1;
    else{
        ids[n]=cnx;
        ++n;
        if(n==max) full=true;
    }

    return 0;
}

int Connector::del(int cnx) {
    //delete a connection

    int shift=-1;
    for(int i=0; i<n; ++i){
        if(ids[i]==cnx){
            shift=i;
            --n;
            full=false;
            break;
        }
    }
    if(shift==-1){
        return 1;
    }
    for(int i=shift; i<n; ++i) ids[i]=ids[i+1];
    ids[n]=-1;

    return 0;
}
