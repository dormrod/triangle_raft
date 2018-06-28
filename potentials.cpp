#include "potentials.h"

//##### BASE POTENTIAL MODEL IN CARTESIAN 2D COORDINATES #####
BasePotentialCart2D::BasePotentialCart2D(){
    //default constructor
}

BasePotentialCart2D::BasePotentialCart2D(vector<int> &bondsIn, vector<int> &anglesIn, vector<int> &fixedIn, vector<int> &interxIn) {
    //set up initial coordinate, bond and angle column vectors
    bonds=bondsIn;
    angles=anglesIn;
    fixed=fixedIn;
    interx=interxIn;
    nBonds=bonds.n/2;
    nAngles=angles.n/3;
    nInterx=interx.n/4;
}

void BasePotentialCart2D::calculateForce(col_vector<double> &crds, col_vector<double> &force) {
    //calculate force from all interactions

    //reset forces
    force=0.0;

    //calculate force from bonds
    int x0, x1, y0, y1; //indices
    for(int i=0; i<nBonds; ++i){
        x0=2*bonds[2*i];
        x1=2*bonds[2*i+1];
        y0=x0+1;
        y1=x1+1;
        bondForce(crds[x0],crds[y0],crds[x1],crds[y1],force[x0],force[y0],force[x1],force[y1],i);
    }

    //calculate force from angles
    int x2, y2; //indices
    for(int i=0; i<nAngles; ++i){
        x0=2*angles[3*i];
        x1=2*angles[3*i+1];
        x2=2*angles[3*i+2];
        y0=x0+1;
        y1=x1+1;
        y2=x2+1;
        angleForce(crds[x0],crds[y0],crds[x1],crds[y1],crds[x2],crds[y2],force[x0],force[y0],force[x1],force[y1],force[x2],force[y2],i);
    }

    //kill forces on fixed atoms
    for(int i=0; i<fixed.n; ++i){
        force[2*fixed[i]]=0.0;
        force[2*fixed[i]+1]=0.0;
    }
}

void BasePotentialCart2D::calculateEnergy(col_vector<double> &crds, double &energy) {
    //calculate energy from all interactions

    //reset energy
    energy=0.0;

    //calculate energy from bonds
    int x0, x1, y0, y1; //indices
    for(int i=0; i<nBonds; ++i){
        x0=2*bonds[2*i];
        x1=2*bonds[2*i+1];
        y0=x0+1;
        y1=x1+1;
        bondEnergy(crds[x0],crds[y0],crds[x1],crds[y1],energy,i);
    }

    //calculate energy from angles
    int x2, y2; //indices
    for(int i=0; i<nAngles; ++i){
        x0=2*angles[3*i];
        x1=2*angles[3*i+1];
        x2=2*angles[3*i+2];
        y0=x0+1;
        y1=x1+1;
        y2=x2+1;
        angleEnergy(crds[x0],crds[y0],crds[x1],crds[y1],crds[x2],crds[y2],energy,i);
    }

    //calculate energy from intersections
    int x3, y3; //indices
    for(int i=0; i<nInterx; ++i){
        x0=2*interx[4*i];
        x1=2*interx[4*i+1];
        x2=2*interx[4*i+2];
        x3=2*interx[4*i+3];
        y0=x0+1;
        y1=x1+1;
        y2=x2+1;
        y3=x3+1;
        interxEnergy(crds[x0],crds[y0],crds[x1],crds[y1],crds[x2],crds[y2],crds[x3],crds[y3],energy);
    }

}

//##### HARMONIC CARTESIAN 2D //
HC2::HC2(){
    //default constructor
}

HC2::HC2(vector<int> &bondsIn, vector<int> &anglesIn, double &bondKIn, double &bondR0In,
         double &angleKIn, double &angleR0In, vector<int> &fixedIn, vector<int> &interxIn):BasePotentialCart2D(bondsIn,anglesIn,fixedIn,interxIn) {
    //construct with single parameter set for bonds and angles
    //turn single values into vectors
    bondK=col_vector<double>(bonds.n);
    bondR0=col_vector<double>(bonds.n);
    angleK=col_vector<double>(angles.n);
    angleR0=col_vector<double>(angles.n);
    bondK=bondKIn;
    bondR0=bondR0In;
    angleK=angleKIn;
    angleR0=angleR0In;
}

HC2::HC2(vector<int> &bondsIn, vector<int> &anglesIn, vector<double> &bondKIn, vector<double> &bondR0In,
         vector<double> &angleKIn, vector<double> &angleR0In, vector<int> &fixedIn, vector<int> &interxIn):BasePotentialCart2D(bondsIn,anglesIn,fixedIn,interxIn) {
    //construct with parameter set for individual bonds and angles
    bondK=bondKIn;
    bondR0=bondR0In;
    angleK=angleKIn;
    angleR0=angleR0In;
}

inline void HC2::bondForce(double &cx0, double &cy0, double &cx1, double &cy1, double &fx0, double &fy0,
                                  double &fx1, double &fy1, int paramRef) {
    //calculate force from single harmonic bond, f=-k(r-r0)
    col_vector<double> f(2);
    f[0]=cx1-cx0;
    f[1]=cy1-cy0;
    double r=sqrt(f[0]*f[0]+f[1]*f[1]);
    double mag=-bondK[paramRef]*(r-bondR0[paramRef])/r;
    f*=mag;
    fx0-=f[0];
    fy0-=f[1];
    fx1+=f[0];
    fy1+=f[1];
}

inline void HC2::angleForce(double &cx0, double &cy0, double &cx1, double &cy1, double &cx2, double &cy2, double &fx0,
                            double &fy0, double &fx1, double &fy1, double &fx2, double &fy2, int paramRef) {
    //calculate force as single harmonic bond, f=-k(r-r0), between outside atoms of angle: neglect central atom
    col_vector<double> f(2);
    f[0]=cx2-cx0;
    f[1]=cy2-cy0;
    double r=sqrt(f[0]*f[0]+f[1]*f[1]);
    double mag=-angleK[paramRef]*(r-angleR0[paramRef])/r;
    f*=mag;
    fx0-=f[0];
    fy0-=f[1];
    fx2+=f[0];
    fy2+=f[1];
}

inline void HC2::bondEnergy(double &cx0, double &cy0, double &cx1, double &cy1, double &e, int paramRef) {
    //calculate energy of a single harmonic bond, U=0.5k(r-r0)^2
    double dx=cx1-cx0;
    double dy=cy1-cy0;
    double r=sqrt(dx*dx+dy*dy);
    e+=0.5*bondK[paramRef]*pow((r-bondR0[paramRef]),2);
}

inline void HC2::angleEnergy(double &cx0, double &cy0, double &cx1, double &cy1, double &cx2, double &cy2, double &e, int paramRef) {
    //calculate energy of single harmonic bond, U=0.5k(r-r0)^2, between outside atoms of angle
    double dx=cx2-cx0;
    double dy=cy2-cy0;
    double r=sqrt(dx*dx+dy*dy);
    e+=0.5*angleK[paramRef]*pow((r-angleR0[paramRef]),2);
}

inline void HC2::interxEnergy(double &cx0, double &cy0, double &cx1, double &cy1, double &cx2, double &cy2, double &cx3,
                              double &cy3, double &e) {
    //calculate energy of single intersection pair, delta function
    bool intersection=properIntersectionLines(cx0,cy0,cx1,cy1,cx2,cy2,cx3,cy3);
    if(intersection) e=numeric_limits<double>::infinity();
}

//##### BASE POTENTIAL MODEL IN CARTESIAN 3D COORDINATES #####
BasePotentialCart3D::BasePotentialCart3D(){
    //default constructor
}

BasePotentialCart3D::BasePotentialCart3D(vector<int> &bondsIn, vector<int> &anglesIn, vector<int> &fixedIn,
                                         vector<int> &interxIn, vector<int> &constrainedIn) {
    //set up initial coordinate, bond and angle column vectors
    bonds=bondsIn;
    angles=anglesIn;
    fixed=fixedIn;
    interx=interxIn;
    constrained=constrainedIn;
    nBonds=bonds.n/2;
    nAngles=angles.n/3;
    nInterx=interx.n/4;
}

void BasePotentialCart3D::calculateForce(col_vector<double> &crds, col_vector<double> &force) {
    //calculate force from all interactions

    //reset forces
    force=0.0;

    //calculate forces from geometry constraints
    int x0, y0, z0; //indices
    for(int i=0; i<constrained.n; ++i){
        x0=3*constrained[i];
        y0=x0+1;
        z0=x0+2;
        constraintForce(crds[x0],crds[y0],crds[z0],force[x0],force[y0],force[z0],i);
    }

    //calculate force from bonds
    int x1, y1, z1; //indices
    for(int i=0; i<nBonds; ++i){
        x0=3*bonds[2*i];
        x1=3*bonds[2*i+1];
        y0=x0+1;
        y1=x1+1;
        z0=x0+2;
        z1=x1+2;
        bondForce(crds[x0],crds[y0],crds[z0],crds[x1],crds[y1],crds[z1],force[x0],force[y0],force[z0],force[x1],force[y1],force[z1],i);
    }

    //calculate force from angles
    int x2, y2, z2; //indices
    for(int i=0; i<nAngles; ++i){
        x0=3*angles[3*i];
        x1=3*angles[3*i+1];
        x2=3*angles[3*i+2];
        y0=x0+1;
        y1=x1+1;
        y2=x2+1;
        z0=x0+2;
        z1=x1+2;
        z2=x2+2;
        angleForce(crds[x0],crds[y0],crds[z0],crds[x1],crds[y1],crds[z1],crds[x2],crds[y2],crds[z2],
                   force[x0],force[y0],force[z0],force[x1],force[y1],force[z1],force[x2],force[y2],force[z2],i);
    }

    //kill forces on fixed atoms
    for(int i=0; i<fixed.n; ++i){
        force[3*fixed[i]]=0.0;
        force[3*fixed[i]+1]=0.0;
        force[3*fixed[i]+2]=0.0;
    }
}

void BasePotentialCart3D::calculateEnergy(col_vector<double> &crds, double &energy) {
    //calculate energy from all interactions

    //reset energy
    energy=0.0;

    //calculate energy from geometry constraints
    int x0, y0, z0; //indices
    for(int i=0; i<constrained.n; ++i){
        x0=3*constrained[i];
        y0=x0+1;
        z0=x0+2;
        constraintEnergy(crds[x0],crds[y0],crds[z0],energy,i);
    }

    //calculate energy from bonds
    int x1, y1, z1; //indices
    for(int i=0; i<nBonds; ++i){
        x0=3*bonds[2*i];
        x1=3*bonds[2*i+1];
        y0=x0+1;
        y1=x1+1;
        z0=x0+2;
        z1=x1+2;
        bondEnergy(crds[x0],crds[y0],crds[z0],crds[x1],crds[y1],crds[z1],energy,i);
    }

    //calculate energy from angles
    int x2, y2, z2; //indices
    for(int i=0; i<nAngles; ++i){
        x0=3*angles[3*i];
        x1=3*angles[3*i+1];
        x2=3*angles[3*i+2];
        y0=x0+1;
        y1=x1+1;
        y2=x2+1;
        z0=x0+2;
        z1=x1+2;
        z2=x2+2;
        angleEnergy(crds[x0],crds[y0],crds[z0],crds[x1],crds[y1],crds[z1],crds[x2],crds[y2],crds[z2],energy,i);
    }

    //calculate energy from intersections
    int x3, y3, z3; //indices
    for(int i=0; i<nInterx; ++i){
        x0=3*interx[4*i];
        x1=3*interx[4*i+1];
        x2=3*interx[4*i+2];
        x3=3*interx[4*i+3];
        y0=x0+1;
        y1=x1+1;
        y2=x2+1;
        y3=x3+1;
        z0=x0+2;
        z1=x1+2;
        z2=x2+2;
        z3=x3+2;
        interxEnergy(crds[x0],crds[y0],crds[z0],crds[x1],crds[y1],crds[z1],crds[x2],crds[y2],crds[z2],crds[x3],crds[y3],crds[z3],energy);
    }

}

//##### HARMONIC CARTESIAN 3D SPHERE CONSTRAINED //
HC3S::HC3S(){
    //default constructor
}

HC3S::HC3S(vector<int> &bondsIn, vector<int> &anglesIn, vector<int> &fixedIn, vector<int> &interxIn,
           vector<int> constrIn, double &bondKIn, double &bondR0In, double &angleKIn, double &angleR0In,
           double &constrKIn, double &constrR0In):BasePotentialCart3D(bondsIn,anglesIn,fixedIn,interxIn,constrIn) {
    //construct with single parameter set for bonds and angles and constraints
    //turn single values into vectors
    bondK=col_vector<double>(bonds.n);
    bondR0=col_vector<double>(bonds.n);
    angleK=col_vector<double>(angles.n);
    angleR0=col_vector<double>(angles.n);
    constraintK=col_vector<double>(constrained.n);
    constraintR0=col_vector<double>(constrained.n);
    bondK=bondKIn;
    bondR0=bondR0In;
    angleK=angleKIn;
    angleR0=angleR0In;
    constraintK=constrKIn;
    constraintR0=constrR0In;
}

HC3S::HC3S(vector<int> &bondsIn, vector<int> &anglesIn, vector<int> &fixedIn, vector<int> &interxIn,
           vector<int> constrIn, vector<double> &bondKIn, vector<double> &bondR0In, vector<double> &angleKIn,
           vector<double> &angleR0In, vector<double> &constrKIn, vector<double> &constrR0In):BasePotentialCart3D(bondsIn,anglesIn,fixedIn,interxIn,constrIn) {
    //construct with parameter set for individual bonds and angles and constraints
    bondK=bondKIn;
    bondR0=bondR0In;
    angleK=angleKIn;
    angleR0=angleR0In;
    constraintK=constrKIn;
    constraintR0=constrR0In;
}

inline void HC3S::constraintForce(double &cx0, double &cy0, double &cz0, double &fx0, double &fy0, double &fz0, int paramRef) {
    //calculate force from point to constraining sphere, f=-k(r-r0)
    col_vector<double> f(3);
    f[0]=cx0;
    f[1]=cy0;
    f[2]=cz0;
    double r=sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
    double mag=-constraintK[paramRef]*(r-constraintR0[paramRef])/r;
    f*=mag;
    fx0+=f[0];
    fy0+=f[1];
    fz0+=f[2];
}

inline void HC3S::bondForce(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &fx0,
                            double &fy0, double &fz0, double &fx1, double &fy1, double &fz1, int paramRef) {
    //calculate force from single harmonic bond, f=-k(r-r0)
    col_vector<double> f(3);
    f[0]=cx1-cx0;
    f[1]=cy1-cy0;
    f[2]=cz1-cz0;
    double r=sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
    double mag=-bondK[paramRef]*(r-bondR0[paramRef])/r;
    f*=mag;
    fx0-=f[0];
    fy0-=f[1];
    fz0-=f[2];
    fx1+=f[0];
    fy1+=f[1];
    fz1+=f[2];
}

inline void HC3S::angleForce(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &cx2,
                             double &cy2, double &cz2, double &fx0, double &fy0, double &fz0, double &fx1, double &fy1,
                             double &fz1, double &fx2, double &fy2, double &fz2, int paramRef) {
    //calculate force as single harmonic bond, f=-k(r-r0), between outside atoms of angle: neglect central atom
    col_vector<double> f(3);
    f[0]=cx2-cx0;
    f[1]=cy2-cy0;
    f[2]=cz2-cz0;
    double r=sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
    double mag=-angleK[paramRef]*(r-angleR0[paramRef])/r;
    f*=mag;
    fx0-=f[0];
    fy0-=f[1];
    fz0-=f[2];
    fx2+=f[0];
    fy2+=f[1];
    fz2+=f[2];
}

inline void HC3S::constraintEnergy(double &cx0, double &cy0, double &cz0, double &e, int paramRef) {
    //calculate energy from point to constraining sphere, U=0.5k(r-r0)^2
    double r=sqrt(cx0*cx0+cy0*cy0+cz0*cz0);
    e+=0.5*constraintK[paramRef]*pow((r-constraintR0[paramRef]),2);
}

inline void HC3S::bondEnergy(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &e,
                             int paramRef) {
    //calculate energy of a single harmonic bond, U=0.5k(r-r0)^2
    double dx=cx1-cx0;
    double dy=cy1-cy0;
    double dz=cz1-cz0;
    double r=sqrt(dx*dx+dy*dy+dz*dz);
    e+=0.5*bondK[paramRef]*pow((r-bondR0[paramRef]),2);
}

inline void HC3S::angleEnergy(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1, double &cx2,
                              double &cy2, double &cz2, double &e, int paramRef) {
    //calculate energy of single harmonic bond, U=0.5k(r-r0)^2, between outside atoms of angle
    double dx=cx2-cx0;
    double dy=cy2-cy0;
    double dz=cz2-cz0;
    double r=sqrt(dx*dx+dy*dy+dz*dz);
    e+=0.5*angleK[paramRef]*pow((r-angleR0[paramRef]),2);
}

inline void HC3S::interxEnergy(double &cx0, double &cy0, double &cz0, double &cx1, double &cy1, double &cz1,
                               double &cx2, double &cy2, double &cz2, double &cx3, double &cy3, double &cz3,
                               double &e) {
    //NOT IMPLEMENTED
//    bool intersection=properIntersectionLines(cx0,cy0,cx1,cy1,cx2,cy2,cx3,cy3);
//    if(intersection) e=numeric_limits<double>::infinity();
}

