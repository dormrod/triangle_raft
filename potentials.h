//Contains potential models.
#ifndef MX2_POTENTIALS_H
#define MX2_POTENTIALS_H

#include <iostream>
#include <cmath>
#include <limits>
#include "col_vector.h"
#include "crd.h"

struct BasePotentialCart2D{
    //base struct for potential model in 2D cartesian coordinates

    int nBonds, nAngles, nRep, nInterx; //number of bonds and angles and repulsions and intersections
    col_vector<int> bonds; //list of bonds
    col_vector<int> angles; //list of angles
    col_vector<int> repulsions; //list of repulsions
    col_vector<int> interx; //list of intersections
    col_vector<int> fixed; //list of fixed atoms

    //constructors
    BasePotentialCart2D();
    BasePotentialCart2D(vector<int> &bondsIn, vector<int> &anglesIn, vector<int> &repIn, vector<int> &fixedIn, vector<int> &interxIn);

    //methods
    void calculateForce(col_vector<double> &crds, col_vector<double> &force);
    void calculateEnergy(col_vector<double> &crds, double &energy);

    //virutal methods
    virtual void bondForce(double &cx0, double &cy0, double &cx1, double &cy1,
                           double &fx0, double &fy0, double &fx1, double &fy1, int paramRef)=0;
    virtual void angleForce(double &cx0, double &cy0, double &cx1, double &cy1, double &cx2, double &cy2,
                            double &fx0, double &fy0, double &fx1, double &fy1, double &fx2, double &fy2, int paramRef)=0;
    virtual void repForce(double &cx0, double &cy0, double &cx1, double &cy1,
                           double &fx0, double &fy0, double &fx1, double &fy1, int paramRef)=0;
    virtual void bondEnergy(double &cx0, double &cy0, double &cx1, double &cy1, double &e, int paramRef)=0;
    virtual void angleEnergy(double &cx0, double &cy0, double &cx1, double &cy1, double &cx2, double &cy2, double &e, int paramRef)=0;
    virtual void repEnergy(double &cx0, double &cy0, double &cx1, double &cy1, double &e, int paramRef)=0;
    virtual void interxEnergy(double &cx0, double &cy0, double &cx1, double &cy1, double &cx2, double &cy2, double &cx3, double &cy3, double &e)=0;
};

struct HC2: public BasePotentialCart2D{
    //harmonic,cartesian, 2D

    //potential information
    col_vector<double> bondK, bondR0, angleK, angleR0; //constant and separation minimum

    //constructors
    HC2();
    HC2(vector<int> &bondsIn, vector<int> &anglesIn, vector<int> &repIn, double &bondKIn, double &bondR0In, double &angleKIn, double &angleR0In, vector<int> &fixedIn, vector<int> &interxIn);
    HC2(vector<int> &bondsIn, vector<int> &anglesIn, vector<int> &repIn, vector<double> &bondKIn, vector<double> &bondR0In, vector<double> &angleKIn, vector<double> &angleR0In, vector<int> &fixedIn, vector<int> &interxIn);

    //overrides for virtual methods
    void bondForce(double &cx0, double &cy0, double &cx1, double &cy1,
                   double &fx0, double &fy0, double &fx1, double &fy1, int paramRef) override;
    void angleForce(double &cx0, double &cy0, double &cx1, double &cy1, double &cx2, double &cy2,
                            double &fx0, double &fy0, double &fx1, double &fy1, double &fx2, double &fy2, int paramRef) override;;
    void repForce(double &cx0, double &cy0, double &cx1, double &cy1,
                          double &fx0, double &fy0, double &fx1, double &fy1, int paramRef) override;;
    void bondEnergy(double &cx0, double &cy0, double &cx1, double &cy1, double &e, int paramRef) override;
    void angleEnergy(double &cx0, double &cy0, double &cx1, double &cy1, double &cx2, double &cy2, double &e, int paramRef) override;
    void repEnergy(double &cx0, double &cy0, double &cx1, double &cy1, double &e, int paramRef) override;
    void interxEnergy(double &cx0, double &cy0, double &cx1, double &cy1, double &cx2, double &cy2, double &cx3, double &cy3, double &e) override;
};

struct HLJC2: public BasePotentialCart2D{
    //harmonic,cartesian, 2D

    //potential information
    col_vector<double> bondK, bondR0, repEpsilon, repR02; //constant and separation minimum

    //constructors
    HLJC2();
    HLJC2(vector<int> &bondsIn, vector<int> &anglesIn, vector<int> &repIn, double &bondKIn, double &bondR0In, double &repEpIn, double &repR0In, vector<int> &fixedIn, vector<int> &interxIn);
    HLJC2(vector<int> &bondsIn, vector<int> &anglesIn, vector<int> &repIn, vector<double> &bondKIn, vector<double> &bondR0In, vector<double> &repEpIn, vector<double> &repR0In, vector<int> &fixedIn, vector<int> &interxIn);

    //overrides for virtual methods
    void bondForce(double &cx0, double &cy0, double &cx1, double &cy1,
                   double &fx0, double &fy0, double &fx1, double &fy1, int paramRef) override;
    void angleForce(double &cx0, double &cy0, double &cx1, double &cy1, double &cx2, double &cy2,
                            double &fx0, double &fy0, double &fx1, double &fy1, double &fx2, double &fy2, int paramRef) override;;
    void repForce(double &cx0, double &cy0, double &cx1, double &cy1,
                          double &fx0, double &fy0, double &fx1, double &fy1, int paramRef) override;;
    void bondEnergy(double &cx0, double &cy0, double &cx1, double &cy1, double &e, int paramRef) override;
    void angleEnergy(double &cx0, double &cy0, double &cx1, double &cy1, double &cx2, double &cy2, double &e, int paramRef) override;
    void repEnergy(double &cx0, double &cy0, double &cx1, double &cy1, double &e, int paramRef) override;
    void interxEnergy(double &cx0, double &cy0, double &cx1, double &cy1, double &cx2, double &cy2, double &cx3, double &cy3, double &e) override;
};

struct BasePotentialCart3D{
    //base struct for potential model in 3D cartesian coordinates, with potential for geometrical constraint

    int nBonds, nAngles, nRep, nInterx; //number of bonds and angles, intersections
    col_vector<int> bonds; //list of bonds
    col_vector<int> angles; //list of angles
    col_vector<int> repulsions; //list of repulsions
    col_vector<int> interx; //list of intersections
    col_vector<int> fixed; //list of fixed atoms
    col_vector<int> constrained; //list of constrained atoms

    //constructors
    BasePotentialCart3D();
    BasePotentialCart3D(vector<int> &bondsIn, vector<int> &anglesIn, vector<int> &repIn, vector<int> &fixedIn, vector<int> &interxIn, vector<int> &constrainedIn);

    //methods
    void calculateForce(col_vector<double> &crds, col_vector<double> &force);
    void calculateEnergy(col_vector<double> &crds, double &energy);

    //virutal methods
    virtual void constraintForce(double &cx0, double &cy0, double &cz0,
                                 double &fx0, double &fy0, double &fz0, int paramRef)=0;
    virtual void bondForce(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1,
                           double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, int paramRef)=0;
    virtual void angleForce(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &cx2, double &cy2, double &cz2,
                            double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, double &fx2, double &fy2, double &fz2, int paramRef)=0;
    virtual void repForce(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1,
                           double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, int paramRef)=0;
    virtual void constraintEnergy(double &cx0, double &cy0, double &cz0, double &e, int paramRef)=0;
    virtual void bondEnergy(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &e, int paramRef)=0;
    virtual void angleEnergy(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &cx2, double &cy2, double &cz2, double &e, int paramRef)=0;
    virtual void repEnergy(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &e, int paramRef)=0;
    virtual void interxEnergy(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &cx2, double &cy2, double &cz2, double &cx3, double &cy3, double &cz3, double &e)=0;
};

struct HC3S: public BasePotentialCart3D{
    //harmonic, cartesian, 3D, sphere constrained

    //potential information
    col_vector<double> bondK, bondR0, angleK, angleR0, constraintR0, constraintK; //constant and separation minimum

    //constructors
    HC3S();
    HC3S(vector<int> &bondsIn, vector<int> &anglesIn, vector<int> &repIn, vector<int> &fixedIn, vector<int> &interxIn, vector<int> constrIn,
         double &bondKIn, double &bondR0In, double &angleKIn, double &angleR0In, double &constrKIn, double &constrR0In);
    HC3S(vector<int> &bondsIn, vector<int> &anglesIn, vector<int> &repIn, vector<int> &fixedIn, vector<int> &interxIn, vector<int> constrIn,
         vector<double> &bondKIn, vector<double> &bondR0In, vector<double> &angleKIn, vector<double> &angleR0In, vector<double> &constrKIn, vector<double> &constrR0In);

    //overrides for virtual methods
    void constraintForce(double &cx0, double &cy0, double &cz0,
                                 double &fx0, double &fy0, double &fz0, int paramRef) override;
    void bondForce(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1,
                           double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, int paramRef) override;
    void angleForce(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &cx2, double &cy2, double &cz2,
                            double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, double &fx2, double &fy2, double &fz2, int paramRef) override;
    void repForce(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1,
                   double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, int paramRef) override;
    void constraintEnergy(double &cx0, double &cy0, double &cz0, double &e, int paramRef) override;
    void bondEnergy(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &e, int paramRef) override;
    void angleEnergy(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &cx2, double &cy2, double &cz2, double &e, int paramRef) override;
    void repEnergy(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &e, int paramRef) override;
    void interxEnergy(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &cx2, double &cy2, double &cz2, double &cx3, double &cy3, double &cz3, double &e) override;
};

struct HLJC3S: public BasePotentialCart3D{
    //harmonic with LJ repulsions, cartesian, 3D, sphere constrained

    //potential information
    col_vector<double> bondK, bondR0, repEpsilon, repR02, constraintR0, constraintK; //constant and separation minimum

    //constructors
    HLJC3S();
    HLJC3S(vector<int> &bondsIn, vector<int> &anglesIn, vector<int> &repIn, vector<int> &fixedIn, vector<int> &interxIn, vector<int> constrIn,
    double &bondKIn, double &bondR0In, double &repEpIn, double &repR0In, double &constrKIn, double &constrR0In);
    HLJC3S(vector<int> &bondsIn, vector<int> &anglesIn, vector<int> &repIn, vector<int> &fixedIn, vector<int> &interxIn, vector<int> constrIn,
    vector<double> &bondKIn, vector<double> &bondR0In, vector<double> &repEpIn, vector<double> &repR0In, vector<double> &constrKIn, vector<double> &constrR0In);

    //overrides for virtual methods
    void constraintForce(double &cx0, double &cy0, double &cz0,
                         double &fx0, double &fy0, double &fz0, int paramRef) override;
    void bondForce(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1,
                   double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, int paramRef) override;
    void angleForce(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &cx2, double &cy2, double &cz2,
                    double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, double &fx2, double &fy2, double &fz2, int paramRef) override;
    void repForce(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1,
                  double &fx0, double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, int paramRef) override;
    void constraintEnergy(double &cx0, double &cy0, double &cz0, double &e, int paramRef) override;
    void bondEnergy(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &e, int paramRef) override;
    void angleEnergy(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &cx2, double &cy2, double &cz2, double &e, int paramRef) override;
    void repEnergy(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &e, int paramRef) override;
    void interxEnergy(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &cx2, double &cy2, double &cz2, double &cx3, double &cy3, double &cz3, double &e) override;
};

#endif //MX2_POTENTIALS_H
