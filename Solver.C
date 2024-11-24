#include "Solver.h"
#include "conservationform.h"

// #define PARALLEL_FOR for (int _i_dummy = 0; _i_dummy < 1; ++_i_dummy) if (parallel.choice == 1) _Pragma("omp parallel for")


void Solver2D::solve(int choice)
{
    check_folder();
    dt = 0.0;
    this->initfunction();
    int step = 0;
    double t = 0.0;
    saveData(t,step); // 保存初始数据 
    auto io_start = std::chrono::high_resolution_clock::now();
    auto io_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> io_duration;
    double total_io_time;
    auto start = std::chrono::high_resolution_clock::now();
    for (;;){
        if (choice == 1){
            timeStep(t);
        }
        getdt();
        t = t + dt;
        step++;
        if(t>=init.t1){
            if(init.timer == 1){
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double, std::milli> duration = end - start;
                std::cout << "\nBlock run time: " << duration.count() << " millisecond" << std::endl;
                std::cout << "Total I/O time: " << total_io_time << " millisecond" << std::endl;
                std::cout << "real computation time: " << duration.count() - total_io_time << " millisecond" << std::endl;
            }
            break;
        }

        applyBoundaryConditions();

        reconstructStates(0);

        computeFluxes(0);

        updateSolution(0);

        reconstructStates(1);
        
        computeFluxes(1);

        updateSolution(1);

        old_u_update();

        io_start = std::chrono::high_resolution_clock::now();
        saveData(step, t);
        io_end = std::chrono::high_resolution_clock::now();
        io_duration = io_end - io_start;
        total_io_time += io_duration.count();
    } 
    


    

}

void Solver2D::old_u_update(){
    int ny = old_u.size();
    int nx = old_u[0].size();
    PARALLEL_FOR
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            new_u[j][i] = old_u[j][i];
        }
    }
}

void Solver2D::updateSolution(int choice) 
{
    int ny = old_u.size();
    int nx = old_u[0].size();

    double dx = (init.x1 - init.x0) / init.npoints_x;
    double dy = (init.y1 - init.y0) / init.npoints_y;


    if (choice == 0){
        PARALLEL_FOR
        for (int j = 0; j < ny ; j++) {
            for (int i = 0; i < nx - 4; i++) {
                // x方向通量差
                for (int k = 0; k < 4; k++) {
                    old_u[j][i+2][k]= old_u[j][i+2][k] -(dt / dx) * (flux_x[j][i+1][k] - flux_x[j][i][k]);
                }
            }
        }
    }
    

    
    if(choice == 1){
        PARALLEL_FOR
        for (int j=0;j<ny-4;j++){
            for (int i=0;i<nx;i++){
                        // y方向通量差
                    for (int k = 0; k < 4; k++) {
                        old_u[j+2][i][k]= old_u[j+2][i][k] -(dt / dy) * (flux_y[j+1][i][k] - flux_y[j][i][k]);
                    }
            }
        }
    }  
}

void Solver2D::computeFluxes(int choice) 
{
    int ny = old_u.size();
    int nx = old_u[0].size();

    flux_x.resize(ny, std::vector<Vector4>(nx-2));
    flux_y.resize(ny-2, std::vector<Vector4>(nx));

    const double gamma = 1.4;
     // 计算x方向的通量
    if(choice ==0){
        PARALLEL_FOR
        for (int j = 0; j < int(rightStates.size()); j++) {
            for (int i = 1; i < int(rightStates[0].size()); i++) {
                Vector4 UL = leftStates[j][i]; // 左侧状态
                Vector4 ULR = rightStates[j][i];      // 左侧状态
                Vector4 UR = rightStates[j][i - 1];      // 右侧状态
                Vector4 URL = leftStates[j][i - 1];      // 右侧状态
                test_x = j;
                test_y = i;
                flux_x[j][i-1] = computeNumericalFlux(UL,ULR, UR,URL, gamma, 0); // 0表示x方向
                
            }
        }
    }
    
    

    // 计算y方向的通量
    if (choice == 1){
        PARALLEL_FOR
        for (int j = 1; j < int(rightStatesy.size()); j++) {
            for (int i = 0; i < int(rightStatesy[0].size()); i++) {
                Vector4 UL = leftStatesy[j][i]; // 下侧状态
                Vector4 ULR = rightStatesy[j][i];      // 上侧状态
                Vector4 UR =  rightStatesy[j - 1][i];      // 上侧状态
                Vector4 URL = leftStatesy[j - 1][i];      // 下侧状态
                test_x = j;
                test_y = i;
                flux_y[j-1][i] = computeNumericalFlux(UL, ULR, UR, URL, gamma, 1); // 1表示y方向
            }

        }
    }
}

Vector4 Solver2D::computeFluxVector(const ConservationForm& U, int direction) 
{
    
    const double gamma = 1.4;

    double rho = U.density;
    if (rho <= 0.0) {
    std::cerr << "Error: Negative or zero density detected at computeFluxVector()." <<rho<< std::endl;
    // std::cout<<init.nonphysical<<std::endl;
    if (init.nonphysical == 1){
            exit(1);
        }
    }
    double mom_x = U.momentum_x;
    double mom_y = U.momentum_y;
    double E = U.energy;

    double u = mom_x / rho;
    double v = mom_y / rho;

    double kinetic_energy = 0.5 * rho * (u * u + v * v);
    double internal_energy = E - kinetic_energy;

    double p = (gamma - 1.0) * internal_energy;
    // std::cout<<"p: "<<p<<std::endl;

    Vector4 flux;

    if (direction == 0) {
        // x方向通量
        flux[0] = mom_x;
        flux[1] = mom_x * u + p;
        flux[2] = mom_x * v;
        flux[3] = (E + p) * u;
    } else {
        // y方向通量
        flux[0] = mom_y;
        flux[1] = mom_x * v;
        flux[2] = mom_y * v + p;
        flux[3] = (E + p) * v;
    }

    return flux;
}

Vector4 Solver2D::computeNumericalFlux(const Vector4& UL,const Vector4& ULR, const Vector4& UR,const Vector4& URL, double gamma, int direction) 
{
    ConservationForm U_L(0.0,0.0,0.0,0.0), U_R(0.0,0.0,0.0,0.0), U_LR(0.0,0.0,0.0,0.0), U_RL(0.0,0.0,0.0,0.0);
    for (int k = 0; k < 4; k++) {
        U_L[k] = UL[k];
        U_R[k] = UR[k];
        U_LR[k] = ULR[k];
        U_RL[k] = URL[k];
    }
    if (U_L.density<=0.0 || U_R.density<=0.0 || U_LR.density<=0.0 || U_RL.density<=0.0){
        std::cerr << "Error: Negative or zero density detected at 163." << std::endl;
        if (init.nonphysical == 1){
            exit(1);
        }
    }
    Vector4 FL = computeFluxVector(U_L, direction);
    Vector4 FR = computeFluxVector(U_R, direction);
    Vector4 FLR = computeFluxVector(U_LR, direction);
    Vector4 FRL = computeFluxVector(U_RL, direction);

    
    Vector4 UL_05;
    Vector4 UR_05;
    for (int k = 0; k < 4; k++) {
        if(direction ==0){
            UL_05[k] = UL[k] - 0.5*(dt/dx)*(FLR[k] - FL[k]);
            UR_05[k] = UR[k] - 0.5*(dt/dx)*(FR[k] - FRL[k]);
            if ((UL[k] - 0.5*(dt/dx)*(FLR[k] - FL[k]))<=0.00001 && k==0){
                std::cout<<"sss:"<<UL[k]<<" "<<FLR[k]<<" "<<FL[k]<<0.5*(dt/dx)*(FLR[k] - FL[k])<<std::endl;
                std::cerr << "Error: Negative or zero density detected at 180." << std::endl;
                if (init.nonphysical == 1){
                    exit(1);
                }
            }
            if((UR[k] - 0.5*(dt/dx)*(FR[k] - FRL[k]))<=0.000001 && k==0){
                std::cerr << "Error: Negative or zero density detected at 184." << std::endl;
                if (init.nonphysical == 1){
                    exit(1);
            }
            }

        }else{
            UL_05[k] = UL[k] - 0.5*(dt/dy)*(FLR[k] - FL[k]);
            UR_05[k] = UR[k] - 0.5*(dt/dy)*(FR[k] - FRL[k]);
        }
        
    }
    ConservationForm temp1(0.0,0.0,0.0,0.0),temp2(0.0,0.0,0.0,0.0);
    for (int k = 0; k < 4; k++) {
        temp1[k] = UL_05[k];
        temp2[k] = UR_05[k];
    }
    if (temp1.density<=0.0 || temp2.density<=0.0){
        std::cerr << "Error: Negative or zero density detected at 187." <<temp1.density<<" "<<test_x<<" "<<test_y<< std::endl;
        if (init.nonphysical == 1){
            exit(1);
        }
    }
    Vector4 FL_05 = computeFluxVector(temp1,direction);
    Vector4 FR_05 = computeFluxVector(temp2,direction);
    ConservationForm U05_05(0.0,0.0,0.0,0.0);
    if(direction == 0){
        for (int k=0;k<4;k++){
            U05_05[k] = 0.5*(UR_05[k] + UL_05[k]) - 0.5*(dt/dx)*(FL_05[k] - FR_05[k]);
        }
    }else{
       for (int k=0;k<4;k++){
            U05_05[k] = 0.5*(UR_05[k] + UL_05[k]) - 0.5*(dt/dy)*(FL_05[k] - FR_05[k]);
        }
    }
    
    if (U05_05.density<=0.0){
        std::cerr << "Error: Negative or zero density detected at 204." << std::endl;
        if (init.nonphysical == 1){
            exit(1);
        }
    }
    Vector4 RI_flux = computeFluxVector(U05_05,direction);
    Vector4 flux;
    for (int k = 0; k < 4; k++) {
            double LF;
        if (direction == 0){
            LF = 0.5*(FR_05[k] + FL_05[k]) + 0.5*(dx/dt)*(UR_05[k] - UL_05[k]);
        }
        else{
            LF = 0.5*(FR_05[k] + FL_05[k]) + 0.5*(dy/dt)*(UR_05[k] - UL_05[k]);
        }
        flux[k] = 0.5*(RI_flux[k] + LF);
        
    }
    return flux;
    
}

void Solver2D::reconstructStates(int choice) 
{

    int ny = old_u.size();
    int nx = old_u[0].size();

    if(choice == 0){
        leftStates.resize(ny, std::vector<Vector4>(nx-2));
        rightStates.resize(ny, std::vector<Vector4>(nx-2));

        PARALLEL_FOR
        for (int j = 0; j < ny; j++) {
            for (int i = 1; i <= nx - 2; i++) {
                Vector4 deltaL, deltaR, slope;
                for (int k = 0; k < 4; k++) {
                    deltaL[k] = old_u[j][i][k] - old_u[j][i - 1][k];
                    deltaR[k] = old_u[j][i + 1][k] - old_u[j][i][k];
                }
                
                for (int k = 0; k < 4; k++) {
                    slope = applyLimiter(deltaL, deltaR,1);
                    leftStates[j][i-1][k] = old_u[j][i][k] - 0.5 * slope[k];
                }
                
                for (int k = 0; k < 4; k++) {
                    slope = applyLimiter(deltaL, deltaR,0);
                    rightStates[j][i-1][k] = old_u[j][i][k] + 0.5 * slope[k];
                    
                }
            }
        }
    }

    if (choice == 1){
        leftStatesy.resize(ny-2, std::vector<Vector4>(nx));
        rightStatesy.resize(ny-2, std::vector<Vector4>(nx));
        PARALLEL_FOR
        for (int j = 1; j <= ny - 2; j++) {
            for (int i = 0; i < nx; i++) {
                
                Vector4 deltaL, deltaR, slope;
                for (int k = 0; k < 4; k++) {
                    deltaL[k] = old_u[j][i][k] - old_u[j - 1][i][k];
                    deltaR[k] = old_u[j + 1][i][k] - old_u[j][i][k];
                }
                

                
                for (int k = 0; k < 4; k++) {
                    slope = applyLimiter(deltaL, deltaR,1);
                    leftStatesy[j-1][i][k] = old_u[j][i][k] - 0.5 * slope[k];
                    
                }
                
                for (int k = 0; k < 4; k++) {
                    slope = applyLimiter(deltaL, deltaR,0);
                    rightStatesy[j-1][i][k] = old_u[j][i][k] + 0.5 * slope[k];
                    
                }
            }
        }
    }
    


    
}

Vector4 Solver2D::applyLimiter(const Vector4& deltaL, const Vector4& deltaR,int direction) 
{
    Vector4 limitedSlope;
    for (int k = 0; k < 4; k++) {
        double a = deltaL[k];
        double b = deltaR[k];
        double di = 0.5*(b+a);
        double r;

        // 防止除零错误
        if (std::abs(b) < 1e-10) {
            limitedSlope[k] =0.0;
            continue;
        } else {
            r = a / b;
        }

        double phi = 0.0;

        // if(r<=0.0){
        //     phi = 0.0;
        // }
        // if(r>0.0){
        //     phi = std::min((2*r)/(1+r),2/(1+r));
        // }
        // limitedSlope[k] = ((2*phi)/(1+phi)) * di;



        // if(r<=0.0){
            
        //     phi = 0.0;
        // }
        // if(r>0.0){
        //     phi = std::min((r*(1+r))/(1+r*r),2.0/(1+r));
        // }
        
        // limitedSlope[k] = phi * di;

        if(r<=0.0){
            phi = 0.0;
        }
        if(r>0.0 && r<=1.0){
            phi = r;
        }
        if(r>1){
            if(direction == 0){
            phi = std::min(1.0,2.0/(1.0+r));}
            else{
                phi = std::min(1.0,(2.0*r)/(1.0+r));
            }
        }
        limitedSlope[k] = phi * di;


        // if (r <= 0.0) {
        //     phi = 0.0;}
        // if(r>0.0 && r<=0.5){
        //     phi = 2*r;
        // }
        // if(r>0.5 && r<=1.0){
        //     phi = 1.0;
        // }
        // else{
        //     phi = std::min(2.0/(1-r), std::min(r, 2.0));
        // }

        // limitedSlope[k] = ((2*phi)/(1+phi)) * di;
        

   
    }
    return limitedSlope;
}

void Solver2D::applyBoundaryConditions() 
{
    int ny = old_u.size();      
    int nx = old_u[0].size();   

    // Boundary condition types: 0 - Extrapolation (default), 1 - Periodic, 2 - Reflective, 3 - Solid Wall

    // Left boundary
    PARALLEL_FOR
    for (int j = 0; j < ny; j++) {
        switch (init.BoundaryCondition) {
            case 1: // Periodic boundary
                old_u[j][0] = old_u[j][nx - 3];
                old_u[j][1] = old_u[j][nx - 2];
                break;
            case 2: // Reflective boundary
                old_u[j][0] = old_u[j][2];
                old_u[j][1] = old_u[j][2];
                old_u[j][0].momentum_x *= -1; // Reverse momentum in x direction
                old_u[j][1].momentum_x *= -1;
                break;
            case 3: // Solid wall boundary
                old_u[j][0] = old_u[j][2];
                old_u[j][1] = old_u[j][2];
                old_u[j][0].momentum_x = 0; // Set velocity in x direction to zero
                old_u[j][0].momentum_y = 0; // Set velocity in y direction to zero
                old_u[j][1].momentum_x = 0;
                old_u[j][1].momentum_y = 0;
                break;
            default: // Default extrapolation boundary
                old_u[j][0] = old_u[j][2];
                old_u[j][1] = old_u[j][2];
        }
    }

    // Right boundary
    PARALLEL_FOR
    for (int j = 0; j < ny; j++) {
        switch (init.BoundaryCondition) {
            case 1: // Periodic boundary
                old_u[j][nx - 1] = old_u[j][2];
                old_u[j][nx - 2] = old_u[j][3];
                break;
            case 2: // Reflective boundary
                old_u[j][nx - 1] = old_u[j][nx - 3];
                old_u[j][nx - 2] = old_u[j][nx - 3];
                old_u[j][nx - 1].momentum_x *= -1; // Reverse momentum in x direction
                old_u[j][nx - 2].momentum_x *= -1;
                break;
            case 3: // Solid wall boundary
                old_u[j][nx - 1] = old_u[j][nx - 3];
                old_u[j][nx - 2] = old_u[j][nx - 3];
                old_u[j][nx - 1].momentum_x = 0; // Set velocity in x direction to zero
                old_u[j][nx - 1].momentum_y = 0; // Set velocity in y direction to zero
                old_u[j][nx - 2].momentum_x = 0;
                old_u[j][nx - 2].momentum_y = 0;
                break;
            default: // Default extrapolation boundary
                old_u[j][nx - 1] = old_u[j][nx - 3];
                old_u[j][nx - 2] = old_u[j][nx - 3];
        }
    }

    // Top boundary
    PARALLEL_FOR
    for (int i = 0; i < nx; i++) {
        switch (init.BoundaryCondition) {
            case 1: // Periodic boundary
                old_u[0][i] = old_u[ny - 3][i];
                old_u[1][i] = old_u[ny - 2][i];
                break;
            case 2: // Reflective boundary
                old_u[0][i] = old_u[2][i];
                old_u[1][i] = old_u[2][i];
                old_u[0][i].momentum_y *= -1; // Reverse momentum in y direction
                old_u[1][i].momentum_y *= -1;
                break;
            case 3: // Solid wall boundary
                old_u[0][i] = old_u[2][i];
                old_u[1][i] = old_u[2][i];
                old_u[0][i].momentum_x = 0; // Set velocity in x direction to zero
                old_u[0][i].momentum_y = 0; // Set velocity in y direction to zero
                old_u[1][i].momentum_x = 0;
                old_u[1][i].momentum_y = 0;
                break;
            default: // Default extrapolation boundary
                old_u[0][i] = old_u[2][i];
                old_u[1][i] = old_u[2][i];
        }
    }

    // Bottom boundary
    PARALLEL_FOR
    for (int i = 0; i < nx; i++) {
        switch (init.BoundaryCondition) {
            case 1: // Periodic boundary
                old_u[ny - 1][i] = old_u[2][i];
                old_u[ny - 2][i] = old_u[3][i];
                break;
            case 2: // Reflective boundary
                old_u[ny - 1][i] = old_u[ny - 3][i];
                old_u[ny - 2][i] = old_u[ny - 3][i];
                old_u[ny - 1][i].momentum_y *= -1; // Reverse momentum in y direction
                old_u[ny - 2][i].momentum_y *= -1;
                break;
            case 3: // Solid wall boundary
                old_u[ny - 1][i] = old_u[ny - 3][i];
                old_u[ny - 2][i] = old_u[ny - 3][i];
                old_u[ny - 1][i].momentum_x = 0; // Set velocity in x direction to zero
                old_u[ny - 1][i].momentum_y = 0; // Set velocity in y direction to zero
                old_u[ny - 2][i].momentum_x = 0;
                old_u[ny - 2][i].momentum_y = 0;
                break;
            default: // Default extrapolation boundary
                old_u[ny - 1][i] = old_u[ny - 3][i];
                old_u[ny - 2][i] = old_u[ny - 3][i];
        }
    }
}

double Solver2D::getdt()
{
    double amax = this->geta();
    double dx = (this->init.x1 - this->init.x0) / this->init.npoints_x;
    double dy = (this->init.y1 - this->init.y0) / this->init.npoints_y;
    double min_dx_dy = std::min(dx, dy);
    dt = (init.C * min_dx_dy) / amax;
    return dt;
}

double Solver2D::geta()
{
    double max_speed = 0.0;
    const double gamma = 1.4;
    int ny = old_u.size();
    int nx = old_u[0].size();
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double rho = old_u[j][i].density;
            double mom_x = old_u[j][i].momentum_x;
            double mom_y = old_u[j][i].momentum_y;
            double E = old_u[j][i].energy;
            // 计算速度分量
            for (int i=0;i<2;i++){
                double V = 0.0;
                if (i == 0){
                    V = mom_x / rho;
                }
                if(i==1){
                    V = mom_y / rho;
                }
                double kinetic_energy = 0.5 * rho * (V*V);
                double internal_energy = E - kinetic_energy;
                double p = (gamma - 1.0) * internal_energy;
                double c = sqrt(gamma * p / rho);
                double total_speed = V + c;
                
                    if (total_speed > max_speed) {
                    max_speed = total_speed;
                }
                
            }
        }
    }
    return max_speed;
}

void Solver2D::initfunction()
{
    int nx = init.npoints_x;
    int ny = init.npoints_y;
    dx = (init.x1 - init.x0) / nx;
    dy = (init.y1 - init.y0) / ny;
    // 调整解向量的大小以包含幽灵单元
    old_u.resize(ny + 4, std::vector<ConservationForm>(nx + 4));
    new_u.resize(ny + 4, std::vector<ConservationForm>(nx + 4));
    // 初始化网格并设置初始条件
    for (int j = 0; j < ny + 4; j++) { 
        double y = init.y0 + (j - 2) * dy; 
        for (int i = 0; i < nx + 4; i++) { 
            double x = init.x0 + (i - 2) * dx; 
            old_u[j][i] = function(x, y);
        }
    }
}

ConservationForm Solver2D::function(double x,double y)
{
    
    double radius = 0.4;
    double distance = sqrt((x - 1)*(x - 1) + (y - 1)*(y - 1));
    double rho, u, v, p;
    if (distance <= radius) {
        rho = 1.0;
        u = 0.0;
        v = 0.0;
        p = 1.0;
    } else {
        rho = 0.125;
        u = 0.0;
        v = 0.0;
        p = 0.1;
    }
    // if(x>0.5 && y>0.5){
    //     rho = 0.5313;
    //     u = 0.0;
    //     v = 0.0;
    //     p = 0.4;
    // }
    // if (x<0.5 && y>0.5){
    //     rho = 1.0;
    //     u = 0.7276;
    //     v = 0.0;
    //     p = 1.0;
    // }
    // if (x<0.5 && y<0.5){
    //     rho = 0.8;
    //     u = 0.0;
    //     v = 0.0;
    //     p = 1.0;
    // }
    // if(x>0.5 && y<0.5){
    //     rho = 1.0;
    //     u = 0.0;
    //     v = 0.7276;
    //     p = 1.0;
    // }

    
    ConservationForm temp(rho,u,v,p);
    return temp;
}

void Solver2D::check_folder(){
    struct stat info;
    if (stat("data", &info) != 0) {
        system("mkdir data");
    } else if (!(info.st_mode & S_IFDIR)) {
        std::cerr << "'data' exists but is not a directory!" << std::endl;
        exit(1);
    }
    std::string folderPath = "./data";  
    try {
            for (const auto& entry : fs::directory_iterator(folderPath)) {
                fs::remove(entry.path());
            }
            std::cout << "All data files have been cleared" << std::endl;
        } catch (const fs::filesystem_error& e) {
            std::cerr << "文件操作出错: " << e.what() << std::endl;
        }
}

void Solver2D::saveData(int step, double t) 
{
    int ny = old_u.size();
    int nx = old_u[0].size();
    // 调试输出，检查矩阵维度
    // std::cout << "Saving data for step " << step << " with dimensions: ny = " << ny << ", nx = " << nx << std::endl;
    // 构建文件名，例如 data/step_0000.csv
    std::ostringstream filename;
    filename << "data/step_" << std::setw(4) << std::setfill('0') << step << ".csv";
    std::ofstream file(filename.str());

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename.str() << std::endl;
        return;
    }

    // std::cout << "File opened successfully: " << filename.str() << std::endl;

    // 写入时间信息
    file << "# Time: " << t << "\n";

    // 写入数据
    for (int j = 2; j < ny - 2; j++) {
        for (int i = 2; i < nx - 2; i++) {
            double rho = old_u[j][i].density;
            double u = old_u[j][i].momentum_x / rho;
            double v = old_u[j][i].momentum_y / rho;

            
            const double gamma = 1.4;
            double E = old_u[j][i].energy;
            
            
            double kinetic_energy = 0.5 * rho * (u * u + v * v);
            double internal_energy = E - kinetic_energy;
            double p = (gamma - 1.0) * internal_energy;

            
            file << rho << "," << u << "," << v << "," << p;
            if (i < nx - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
    // std::cout << "Data saved successfully to " << filename.str() << std::endl;
}

void Solver2D::timeStep(double t){
    // 打印一行删除一行

    // 计算进度百分比
    double progress = (t / init.t1) * 100.0;

    // 生成进度条
    int barWidth = 50;
    std::cout << "\r进度: [";
    int pos = barWidth * progress / 100;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress) << " %\t" << std::flush;
}

Solver2D::Solver2D(InitState init,Parallel parallel):init(init),parallel(parallel){}