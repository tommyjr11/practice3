#include "conservationform.h"

ConservationForm::ConservationForm(double rho, double v_x, double v_y, double p)
        : density(rho), momentum_x(rho*v_x), momentum_y(rho*v_y), energy(0.5*rho*(v_x*v_x+v_y*v_y) + p/(1.4-1.0)) {}


double& ConservationForm::operator[](int index) {
    if (index == 0) return density;
    else if (index == 1) return momentum_x;
    else if (index == 2) return momentum_y;
    else if (index == 3) return energy;
    else throw std::out_of_range("Index out of range");
}


const double& ConservationForm::operator[](int index) const {
    if (index == 0) return density;
    else if (index == 1) return momentum_x;
    else if (index == 2) return momentum_y;
    else if (index == 3) return energy;
    else throw std::out_of_range("Index out of range");
}