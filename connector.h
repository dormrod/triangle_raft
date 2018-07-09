#ifndef MX2_CONNECTOR_H
#define MX2_CONNECTOR_H

#include <iostream>

using namespace std;

struct Connector {
    //holds connectivity information

    //container variables
    int n, max; //number of connections, maximum number of connections
    int *ids; //list of connections
    bool full; //whether container is full

    //constructors, destructors, overloaded operators
    Connector();
    Connector(int maxN);
    ~Connector();
    Connector(const Connector &source);
    Connector& operator=(const Connector &source);
    friend ostream& operator<<(ostream &output, const Connector &source) {
        for (int i = 0; i < source.max; ++i) output << source.ids[i] << " ";
        return output;
    };

    //methods
    int add(int cnx); //add a connection
    int del(int cnx); //delete a connection
};

#endif //MX2_CONNECTOR_H
