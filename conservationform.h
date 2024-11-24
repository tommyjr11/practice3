#ifndef conservationform_H
#define conservationform_H
#include <vector>
#include <string>
#include <stdexcept>


class ConservationForm {
    public:
        double density;
        double momentum_x;
        double momentum_y;
        double energy;
        ConservationForm(double rho, double v_x, double v_y, double p);
        ConservationForm() : density(0.0), momentum_x(0.0), momentum_y(0.0), energy(0.0) {}
        double& operator[](int index); 
        const double& operator[](int index) const;

};


#endif