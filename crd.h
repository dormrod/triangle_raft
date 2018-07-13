#ifndef MX2_CRD_H
#define MX2_CRD_H

#include <iostream>
#include <cmath>

using namespace std;

struct Cart2D {
    //2 dimensional coordinate

    //x,y coordinates and initialisation
    double x, y;
    Cart2D();
    Cart2D(double xInit, double yInit);

    //overload operators
    Cart2D operator=(const Cart2D &c){
        if(this==&c) return *this;
        this->x=c.x;
        this->y=c.y;
        return *this;
    }
    Cart2D operator+(const Cart2D &c){
        Cart2D crd;
        crd.x=this->x+c.x;
        crd.y=this->y+c.y;
        return crd;
    }
    Cart2D operator-(const Cart2D &c){
        Cart2D crd;
        crd.x=this->x-c.x;
        crd.y=this->y-c.y;
        return crd;
    }
    Cart2D operator*(const double &k){
        Cart2D crd;
        crd.x=this->x*k;
        crd.y=this->y*k;
        return crd;
    }
    double operator*(const Cart2D &c){
        double dotProduct=0.0;
        dotProduct+=this->x*c.x;
        dotProduct+=this->y*c.y;
        return dotProduct;
    }
    Cart2D operator/(const double &k){
        Cart2D crd;
        crd.x=this->x/k;
        crd.y=this->y/k;
        return crd;
    }
    Cart2D operator+=(const Cart2D &c){
        this->x=this->x+c.x;
        this->y=this->y+c.y;
        return *this;
    }
    Cart2D operator-=(const Cart2D &c){
        this->x=this->x-c.x;
        this->y=this->y-c.y;
        return *this;
    }
    Cart2D operator*=(const double &k){
        this->x=this->x*k;
        this->y=this->y*k;
        return *this;
    }
    Cart2D operator/=(const double &k){
        this->x=this->x/k;
        this->y=this->y/k;
        return *this;
    }

    //additional methods
    double norm();
    double normSq();
    void normalise();
    void normalise(double &norm);
};

struct Cart3D {
    //3 dimensional coordinate

    //x,y,z coordinates and initialisation
    double x, y, z;
    Cart3D();
    Cart3D(double xInit, double yInit, double zInit=0.0);
    Cart3D(Cart2D c2d, double zInit=0.0);

    //overload operators
    Cart3D operator=(const Cart3D &c){
        if(this==&c) return *this;
        this->x=c.x;
        this->y=c.y;
        this->z=c.z;
        return *this;
    }
    Cart3D operator=(const Cart2D &c){
        this->x=c.x;
        this->y=c.y;
        this->z=0.0;
        return *this;
    }
    Cart3D operator+(const Cart3D &c){
        Cart3D crd;
        crd.x=this->x+c.x;
        crd.y=this->y+c.y;
        crd.z=this->z+c.z;
        return crd;
    }
    Cart3D operator-(const Cart3D &c){
        Cart3D crd;
        crd.x=this->x-c.x;
        crd.y=this->y-c.y;
        crd.z=this->z-c.z;
        return crd;
    }
    Cart3D operator*(const double &k){
        Cart3D crd;
        crd.x=this->x*k;
        crd.y=this->y*k;
        crd.z=this->z*k;
        return crd;
    }
    double operator*(const Cart3D &c){
        double dotProduct=0.0;
        dotProduct+=this->x*c.x;
        dotProduct+=this->y*c.y;
        dotProduct+=this->z*c.z;
        return dotProduct;
    }
    Cart3D operator/(const double &k){
        Cart3D crd;
        crd.x=this->x/k;
        crd.y=this->y/k;
        crd.z=this->z/k;
        return crd;
    }
    Cart3D operator+=(const Cart3D &c){
        this->x=this->x+c.x;
        this->y=this->y+c.y;
        this->z=this->z+c.z;
        return *this;
    }
    Cart3D operator-=(const Cart3D &c){
        this->x=this->x-c.x;
        this->y=this->y-c.y;
        this->z=this->z-c.z;
        return *this;
    }
    Cart3D operator*=(const double &k){
        this->x=this->x*k;
        this->y=this->y*k;
        this->z=this->z*k;
        return *this;
    }
    Cart3D operator/=(const double &k){
        this->x=this->x/k;
        this->y=this->y/k;
        this->z=this->z/k;
        return *this;
    }

    //additional methods
    double norm();
    double normSq();
    void normalise();
    void normalise(double &norm);
    Cart2D xyProjection();
};

//##### LINE INTERSECTIONS - IMPLEMENTED USING COMPUTATIONAL GEOMETRY IN C #######
inline double signedAreaSqTriangle(double &x0, double &y0, double &x1, double &y1, double &x2, double &y2){
    //returns square of the signed area of a triangle defined by three points
    return x0*y1-y0*x1+y0*x2-x0*y2+x1*y2-x2*y1;
}
inline bool leftTriangle(double &x0, double &y0, double &x1, double &y1, double &x2, double &y2){
    //returns true if the signed area squared of a triangle is positive
    return signedAreaSqTriangle(x0,y0,x1,y1,x2,y2)>0.0;
}
inline bool collinearPoints(double &x0, double &y0, double &x1, double &y1, double &x2, double &y2, double threshold=1e-6){
    //returns true if points are collinear within given threshold
    return fabs(signedAreaSqTriangle(x0,y0,x1,y1,x2,y2))<threshold;
}
inline bool properIntersectionLines(double &x0, double &y0, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3){
    //returns true if lines properly intersect i.e. no points collinear
    if(collinearPoints(x0,y0,x1,y1,x2,y2) || collinearPoints(x0,y0,x1,y1,x3,y3) || collinearPoints(x2,y2,x3,y3,x0,y0) || collinearPoints(x2,y2,x3,y3,x1,y1)) return false;
    return (leftTriangle(x0,y0,x1,y1,x2,y2)^leftTriangle(x0,y0,x1,y1,x3,y3)) && (leftTriangle(x2,y2,x3,y3,x0,y0)^leftTriangle(x2,y2,x3,y3,x1,y1));
}
inline bool between(double &x0, double &y0, double &x1, double &y1, double &x2, double &y2){
    //checks betweeness of points
    if(!collinearPoints(x0,y0,x1,y1,x2,y2)) return false;
    if(x0!=x1) return ((x0 <= x2) && (x2 <= x1)) || ((x0 >= x2) && (x2 >= x1));
    else return ((y0 <= y2) && (y2 <= y1)) || ((y0 >= y2) && (y2 >= y1));
}
inline bool lineIntersect2D(double &x0, double &y0, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3){
    //returns true if proper or improper intersection
    if(properIntersectionLines(x0,y0,x1,y1,x2,y2,x3,y3)) return true;
    else if(between(x0,y0,x1,y1,x2,y2) || between(x0,y0,x1,y1,x3,y3) || between(x2,y2,x3,y3,x0,y0) || between(x2,y2,x3,y3,x1,y1)) return true;
    return false;
}

#endif //MX2_CRD_H
