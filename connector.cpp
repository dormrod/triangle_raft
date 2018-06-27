#include "connector.h"

//##### CONNECTOR #####
Connector::Connector() {
    //default constructor
    max=0;
}

Connector::Connector(int maxN){
    //constructor
    n=0;
    full=false;
    max=maxN;
    cnxs=new int[max]();
}

Connector::~Connector() {
    //destructor
    delete[] cnxs;
}

Connector::Connector(const Connector &source) {
    //copy constructor

    //shallow copies
    n=source.n;
    full=source.full;
    max=source.max;

    //deep copies
    cnxs=new int[max]();
    for(int i=0; i<n; ++i) cnxs[i]=source.cnxs[i];
}

Connector& Connector::operator=(const Connector &source) {
    //overload assignment operator

    if (this == &source) return *this;

    //shallow copies
    n=source.n;
    full=source.full;
    max=source.max;

    //deep copies
    cnxs=new int[max]();
    for(int i=0; i<n; ++i) cnxs[i]=source.cnxs[i];

    return *this;
}

int Connector::add(int cnx) {
    //add a connection if not full, check if full after addition

    if(full) return 1;
    else{
        cnxs[n]=cnx;
        ++n;
        if(n==max) full=true;
    }

    return 0;
}