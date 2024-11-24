#ifndef SOLVER_H
#define SOLVER_H
#include "Solver.h"
#include "conservationform.h"
#include "initstate.h"
#include "Parallel.h"
#include <sstream>   
#include <iomanip>  
#include <fstream>   
#include <sys/stat.h> 
#include <cstdlib>   
#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>
#include <chrono> 
#include <omp.h>
#include <array>
#include <cmath>


typedef std::array<double, 4> Vector4;
namespace fs = std::filesystem;
class Solver2D{
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
        int test_x;
        int test_y;
        // std::vector<Vector4> r_flux;
        Solver2D(InitState init,Parallel parallel);
        ConservationForm function(double x,double y);

        // 斜率限制器
        Vector4 applyLimiter(const Vector4& deltaL, const Vector4& deltaR,int direction);
        // 重构左右状态
        void reconstructStates(int choice);
        // 计算数值通量
        Vector4 computeNumericalFlux(const Vector4& UL,const Vector4& ULR, const Vector4& UR,const Vector4& URL, double gamma, int direction);
        // 计算通量
        void computeFluxes(int choice);

        void old_u_update();

        Vector4 computeFluxVector(const ConservationForm& U, int direction);
        void updateSolution(int choice);
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