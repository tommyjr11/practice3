#include "Solver_HLLC.h"
#include "conservationform.h"

void Solver_HLLC_2D::solve(int choice){
    check_folder();
    dt = 0.0;
    initfunction();
    int step = 0;
    double t = 0.0;
    saveData(t,step);
    auto start = std::chrono::high_resolution_clock::now();
    auto io_start = std::chrono::high_resolution_clock::now();
    auto io_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> io_duration;
    double total_io_time;
    for(;;){
        if (choice == 1){
            timeStep(t);
        }
        getdt();
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
        t += dt;
        step++;
        
        applyBoundaryConditions();
        
        updateSolution();
        
        io_start = std::chrono::high_resolution_clock::now();
        saveData(step, t);
        io_end = std::chrono::high_resolution_clock::now();
        io_duration = io_end - io_start;
        total_io_time += io_duration.count();
    }




}
Vector4 Solver_HLLC_2D::computeFlux(const ConservationForm& U, char direction){
    double gamma = 1.4;
    double rho = U.density;
    if(rho <= 0.0){
        std::cerr << "Error: Negative or zero density detected at 37." << std::endl;
        if (init.nonphysical == 1){
            exit(1);
        }
    }
    double u = U.momentum_x / rho;
    double v = U.momentum_y / rho;
    double E = U.energy;
    double p = (gamma -1.0)*(E - 0.5*rho*(u*u + v*v));

    Vector4 flux;
    if(direction == 'x'){
        flux[0] = rho * u;
        flux[1] = rho * u * u + p;
        flux[2] = rho * u * v;
        flux[3] = (E + p) * u;
    }
    else{ 
        flux[0] = rho * v;
        flux[1] = rho * u * v;
        flux[2] = rho * v * v + p;
        flux[3] = (E + p) * v;
    }
    return flux;
}

Vector4 Solver_HLLC_2D::computeHLLCFlux(const ConservationForm& UL, const ConservationForm& UR, char direction){
    // 计算左状态
    double gamma = 1.4;
    double rhoL = UL.density;
    double uL = UL.momentum_x / rhoL;
    double vL = UL.momentum_y / rhoL;
    double EL = UL.energy;
    double pL = (gamma -1.0)*(EL - 0.5*rhoL*(uL*uL + vL*vL));
    
    // 计算右状态
    double rhoR = UR.density;
    double uR = UR.momentum_x / rhoR;
    double vR = UR.momentum_y / rhoR;
    double ER = UR.energy;
    double pR = (gamma -1.0)*(ER - 0.5*rhoR*(uR*uR + vR*vR));

    // 计算声速
    double cL = sqrt(gamma * pL / rhoL);
    double cR = sqrt(gamma * pR / rhoR);

    // 估计波速
    double SL;
    double SR;
    if(direction == 'x'){
        SL = std::min(uL - cL, uR - cR);
        SR = std::max(uL + cL, uR + cR);
    }
    else{
        SL = std::min(vL - cL, vR - cR);
        SR = std::max(vL + cL, vR + cR);
    }
    
    
    
    double Sstar;
    if(direction == 'x'){
        Sstar = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR))/(rhoL * (SL - uL) - rhoR * (SR - uR));
    }
    else{
        Sstar = (pR - pL + rhoL * vL * (SL - vL) - rhoR * vR * (SR - vR))/(rhoL * (SL - vL) - rhoR * (SR - vR));
    }
    
    Vector4 UhllcL;
    Vector4 UhllcR;
    if (direction == 'x'){
        UhllcL[0] = rhoL * ((SL - uL) / (SL - Sstar));
        UhllcL[1] = rhoL * ((SL - uL) / (SL - Sstar))*Sstar;
        UhllcL[2] = rhoL * ((SL - uL) / (SL - Sstar))*vL;
        UhllcL[3] = (rhoL * ((SL - uL) / (SL - Sstar)))*((EL / rhoL) + (Sstar-uL)*(Sstar + pL / (rhoL * (SL - uL))));
        UhllcR[0] = rhoR * ((SR - uR) / (SR - Sstar));
        UhllcR[1] = rhoR * ((SR - uR) / (SR - Sstar))*Sstar;
        UhllcR[2] = rhoR * ((SR - uR) / (SR - Sstar))*vR;
        UhllcR[3] = (rhoR * ((SR - uR) / (SR - Sstar)))*((ER / rhoR) + (Sstar-uR)*(Sstar + pR / (rhoR * (SR - uR))));
    }
    else{
        UhllcL[0] = rhoL * ((SL - vL) / (SL - Sstar));
        UhllcL[1] = rhoL * ((SL - vL) / (SL - Sstar))*uL;
        UhllcL[2] = rhoL * ((SL - vL) / (SL - Sstar))*Sstar;
        UhllcL[3] = (rhoL * ((SL - vL) / (SL - Sstar)))*((EL / rhoL) + (Sstar-vL)*(Sstar + pL / (rhoL * (SL - vL))));
        UhllcR[0] = rhoR * ((SR - vR) / (SR - Sstar));
        UhllcR[1] = rhoR * ((SR - vR) / (SR - Sstar))*uR;
        UhllcR[2] = rhoR * ((SR - vR) / (SR - Sstar))*Sstar;
        UhllcR[3] = (rhoR * ((SR - vR) / (SR - Sstar)))*((ER / rhoR) + (Sstar-vR)*(Sstar + pR / (rhoR * (SR - vR))));
    }
    
    
    // 判断是否存在星形状态
    if (SL >= 0 ){
        // 计算通量左边
        
        Vector4 FL = computeFlux(UL, direction);
        
        return FL;
    }
    else if (SL <0 && 0 <= Sstar){
        Vector4 FL = computeFlux(UL, direction);
        Vector4 PartR;
        for (int i = 0; i < 4; i++) {
            PartR[i] = SL * (UhllcL[i] - UL[i]);
        }
        Vector4 flux;
        for (int i = 0; i < 4; i++) {
            flux[i] = FL[i] + PartR[i];
        }
        return flux;
    }
    else if(Sstar < 0 && 0 <= SR){
        Vector4 FR = computeFlux(UR, direction);
        Vector4 PartL;
        for (int i = 0; i < 4; i++) {
            PartL[i] = SR * (UhllcR[i] - UR[i]);
        }
        Vector4 flux;
        for (int i = 0; i < 4; i++) {
            flux[i] = FR[i] + PartL[i];
        }
        return flux;
    }
    else if (0>SR){
        Vector4 FR = computeFlux(UR, direction);
        return FR;
    }
    return Vector4();
}

void Solver_HLLC_2D::updateSolution(){
    int ny = old_u.size();
    int nx = old_u[0].size();
    

    flux_x.resize(ny, std::vector<Vector4>(nx-1));
    flux_y.resize(ny-1, std::vector<Vector4>(nx));
    // 计算x方向通量
    PARALLEL_FOR
    for(int j =0; j < ny; j++){
        for(int i =0; i < nx-1; i++){
            // 左右两个状态
            ConservationForm UL = old_u[j][i];
            ConservationForm UR = old_u[j][i+1];
            // 计算HLLC通量
            Vector4 flux = computeHLLCFlux(UL, UR, 'x');
            flux_x[j][i] = flux;
        }
    }


    // 更新新解 x 方向
    PARALLEL_FOR
    for(int j = 0; j < ny; j++){
        for(int i = 1; i < nx-1; i++){
            for (int k=0;k<4;k++){
                old_u[j][i][k] = old_u[j][i][k] - (dt/dx)*(flux_x[j][i][k] - flux_x[j][i-1][k]);
            }
            // 检查密度和压力
        }
    }



    // 计算y方向通量
    PARALLEL_FOR
    for(int j =0; j < ny-1; j++){
        for(int i =0; i < nx; i++){
            // 上下两个状态
            ConservationForm UL = old_u[j][i];
            ConservationForm UR = old_u[j+1][i];
            // 计算HLLC通量
            Vector4 flux = computeHLLCFlux(UL, UR, 'y');
            flux_y[j][i] = flux;
        }
    }
   
    // 更新新解
    PARALLEL_FOR
    for(int j =1; j < ny-1; j++){
        for(int i =0; i < nx; i++){
            for (int k=0;k<4;k++){
            old_u[j][i][k] = old_u[j][i][k] - (dt/dy)*(flux_y[j][i][k] - flux_y[j-1][i][k]);
            }
        }
    }

    PARALLEL_FOR
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            new_u[j][i] = old_u[j][i];
        }
    }
}

void Solver_HLLC_2D::saveData(int step, double t) 
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
    for (int j = 1; j < ny - 1; j++) {
        for (int i = 1; i < nx - 1; i++) {
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

ConservationForm Solver_HLLC_2D::function(double x,double y){
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

void Solver_HLLC_2D::applyBoundaryConditions() 
{
    int ny = old_u.size();      
    int nx = old_u[0].size();   

    // Boundary condition types: 0 - Extrapolation (default), 1 - Periodic, 2 - Reflective, 3 - Solid Wall

    // Left boundary
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

void Solver_HLLC_2D::timeStep(double t){
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

void Solver_HLLC_2D::initfunction()
{
    int nx = init.npoints_x;
    int ny = init.npoints_y;
    dx = (init.x1 - init.x0) / nx;
    dy = (init.y1 - init.y0) / ny;
    // 调整解向量的大小以包含幽灵单元
    old_u.resize(ny + 2, std::vector<ConservationForm>(nx + 2));
    new_u.resize(ny + 2, std::vector<ConservationForm>(nx + 2));
    // 初始化网格并设置初始条件
    for (int j = 0; j < ny + 2; j++) { 
        double y = init.y0 + (j - 1) * dy; 
        for (int i = 0; i < nx + 2; i++) { 
            double x = init.x0 + (i - 1) * dx; 
            old_u[j][i] = function(x, y);
        }
    }

   
}

double Solver_HLLC_2D::getdt()
{
    double amax = this->geta();
    double dx = (this->init.x1 - this->init.x0) / this->init.npoints_x;
    double dy = (this->init.y1 - this->init.y0) / this->init.npoints_y;
    double min_dx_dy = std::min(dx, dy);
    dt = (init.C * min_dx_dy) / amax;
    return dt;
}

double Solver_HLLC_2D::geta()
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

void Solver_HLLC_2D::check_folder(){
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

Solver_HLLC_2D::Solver_HLLC_2D(InitState init,Parallel parallel)
        : init(init),parallel(parallel) {}