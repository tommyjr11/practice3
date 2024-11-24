#include "initstate.h"

InitState::InitState(double to, double t1, double npoints_x, double npoints_y, double C, double x0, double x1,double y0,double y1,int BoundaryCondition, int nonphysical,int timer)
        : t0(to), t1(t1), npoints_x(npoints_x),npoints_y(npoints_y),C(C), x1(x1), x0(x0),y0(y0),y1(y1),BoundaryCondition(BoundaryCondition),nonphysical(nonphysical),timer(timer){}