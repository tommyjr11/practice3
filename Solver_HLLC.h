#ifndef SOLVER_HLLC_H
#define SOLVER_HLLC_H
#include "conservationform.h"
#include "initstate.h"
#include "Parallel.h"
#include <omp.h>
#include <sstream>   
#include <iomanip>  
#include <fstream>   
#include <sys/stat.h> 
#include <cstdlib>   
#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>
#include <cmath>
#include <array>
typedef std::array<double, 4> Vector4;


namespace fs = std::filesystem;
class Solver_HLLC_2D{
    public:
        InitState init;
        Parallel parallel;
        std::vector<std::vector<ConservationForm>> old_u;
        std::vector<std::vector<ConservationForm>> new_u;
        std::vector<std::vector<std::vector<Vector4>>> gradients;
        std::vector<std::vector<Vector4>> leftStates;
        std::vector<std::vector<Vector4>> rightStates;
        std::vector<std::vector<Vector4>> leftStatesy;
        std::vector<std::vector<Vector4>> rightStatesy;
        // 用于存储界面上的通量
        std::vector<std::vector<Vector4>> flux_x;
        std::vector<std::vector<Vector4>> flux_y;
        double dt;
        double dx;
        double dy;
        std::vector<Vector4> r_flux;



        Solver_HLLC_2D(InitState init,Parallel parallel);
        ConservationForm function(double x,double y);

        
        Vector4 computeHLLCFlux(const ConservationForm& UL, const ConservationForm& UR, char direction);
        Vector4 computeFlux(const ConservationForm& U, char direction);
        void updateSolution();
        void applyBoundaryConditions();
        void saveData(int step, double t);
        void check_folder();

        void initfunction();
        double getdt();
        double geta();
        void solve(int choice);
        void timeStep(double t);



};

#endif