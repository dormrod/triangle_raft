#include "crd.h"

Cart2D::Cart2D() {
    //coordinate at origin
    x=0.0;
    y=0.0;
    return;
}

Cart2D::Cart2D(double xInit, double yInit) {
    //coordinate at specified point
    x=xInit;
    y=yInit;
    return;
}

double Cart2D::norm() {
    //get length of equivalent vector

    return sqrt(x*x+y*y);
}

double Cart2D::normSq() {
    //get squared length of equivalent vector

    return x*x+y*y;
}

void Cart2D::normalise() {
    //make equivalent vector unit

    double rn=1.0/sqrt(x*x+y*y);
    x*=rn;
    y*=rn;
}

void Cart2D::normalise(double &norm) {
    //make equivalent vector unit, pass original length

    norm=sqrt(x*x+y*y);
    double rn=1.0/norm;
    x*=rn;
    y*=rn;
}


Cart3D::Cart3D() {
    //coordinate at origin
    x=0.0;
    y=0.0;
    z=0.0;
    return;
}

Cart3D::Cart3D(double xInit, double yInit, double zInit) {
    //coordinate at specified point
    x=xInit;
    y=yInit;
    z=zInit;
    return;
}

Cart3D::Cart3D(Cart2D c2d, double zInit) {
    //coordinate from 2D coordinate
    x=c2d.x;
    y=c2d.y;
    z=zInit;
}

double Cart3D::norm() {
    //get length of equivalent vector

    return sqrt(x*x+y*y+z*z);
}

double Cart3D::normSq() {
    //get squared length of equivalent vector

    return x*x+y*y+z*z;
}

void Cart3D::normalise() {
    //make equivalent vector unit

    double rn=1.0/sqrt(x*x+y*y+z*z);
    x*=rn;
    y*=rn;
    z*=rn;
}

void Cart3D::normalise(double &norm) {
    //make equivalent vector unit, pass original length

    norm=sqrt(x*x+y*y+z*z);
    double rn=1.0/norm;
    x*=rn;
    y*=rn;
    z*=rn;
}

Cart2D Cart3D::xyProjection() {
    //return x and y component

    return Cart2D(x,y);
}