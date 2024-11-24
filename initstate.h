#ifndef INITSTATE_H
#define INITSTATE_H

#include <vector>
#include <string>
class InitState {
public:
    double t0, t1, npoints_x, npoints_y, C, x1, x0,y0,y1;
    int BoundaryCondition;
    int nonphysical;
    int timer;
    InitState(double to, double t1, double npoints_x, double npoints_y, double C, double x0, double x1,double y0,double y1,int BoundaryCondition,int nonphysical,int timer);
};

#endif