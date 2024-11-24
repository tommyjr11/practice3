# README

## Introduction
This project implements numerical methods for solving two-dimensional hyperbolic conservation laws. It includes two numerical solvers:
    - **HLLC Solver**  
    - **SLIC Solver**
The code is written in C++ with support for parallel execution using OpenMP. Visualization of the simulation results is provided through a Python script.

## Requirements
### Software Requirements

- **C++ Compiler**: `g++-14` and a compatible C++17 compiler with OpenMP support.
- **Make**: For building the project using the provided Makefile.
- **Python 3**: For visualization of the results.
- **Python Libraries**:
  - `matplotlib`: For plotting.
  - `numpy`: For numerical operations.
  - `scipy` (optional): If additional scientific computations are needed.

### Clone the Repository

```bash
git clone [repository URL]
cd [repository directory]
```

### Build the Project

A shell script `1.sh` is provided for cleaning, building, and running the code. It contains:

```bash
#!/bin/bash
make clean
make
./bin/1
```
This will compile the source files and place the executable in the `bin` directory.

### Makefile Details

The provided `Makefile` uses `g++-14` and includes flags for C++17 standard, all warnings enabled, and OpenMP support:

```makefile
CXX = g++-14
CXXFLAGS = -std=c++17 -Wall -fopenmp
```

Source files are:

- `conservationform.C`
- `initstate.C`
- `Solver.C`
- `Solver_HLLC.C`
- `Parallel.C`
- `main.C`

The executable is generated at `./bin/1`.

## Running the Code

After building the project, you can run the simulation using:

```bash
./bin/1
```

Alternatively, use the provided script:

```bash
./1.sh
```

### Program Parameters

The simulation parameters are defined within the code and can be adjusted as needed. Key parameters include:

1. **Initial Time**: Start time of the simulation.
2. **End Time**: Finish time of the simulation.
3. **Grid Points in X-direction (Nx)**: Number of grid points along the x-axis.
4. **Grid Points in Y-direction (Ny)**: Number of grid points along the y-axis.
5. **CFL Number**: Courant–Friedrichs–Lewy (CFL) condition number for numerical stability.
6. **X-direction Start Point**: Minimum value of x in the simulation domain.
7. **X-direction End Point**: Maximum value of x in the simulation domain.
8. **Y-direction Start Point**: Minimum value of y in the simulation domain.
9. **Y-direction End Point**: Maximum value of y in the simulation domain.
10. **Boundary Condition Code**:
    - `0`: Default boundary conditions.
    - `1`: Periodic boundary conditions.
    - `2`: Reflective boundary conditions.
    - `3`: Solid wall boundary conditions.
11. **Non-physical Point Detection**: Enable/disable detection of non-physical values.

These parameters can be modified in the `main.C` file or the relevant configuration files/functions.

### Selecting the Solver

In the `main.C` file, you can choose between the two solvers:

- **HLLC Solver**: Use the `Solver_HLLC_2D` class.
- **SLIC Solver**: Use the `Solver2D` class.

Instantiate the desired solver in the `main` function.

## Customization

### Initial Conditions

The initial conditions for the simulation are defined in the `initstate.C` file or within each solver class in a function named `function`. You can modify this function to set up different initial conditions as per your problem requirements.

### Boundary Conditions

Boundary conditions are set using an integer code (as described above) and can be adjusted in the simulation setup section of the code. Ensure that the selected boundary condition code matches the physical scenario you are modeling.

### Limiter Functions

Limiter functions control numerical diffusion and are essential for capturing shocks and discontinuities without introducing spurious oscillations. You can change the limiter function in the `applylimiter` function within each solver class. Some limiter functions may be commented out due to potential errors. Uncomment and modify them carefully.

### Conservation Form

The `conservationform` class initializes the conservation variables. By default, calling `conservationform()` initializes all variables to `0.0`. If parameters are provided, it defaults to the primitive form and automatically converts to the conservation form. Future versions may allow for selection between different forms.

### Detection of Non-physical Points

The code includes an option to detect non-physical points (e.g., negative densities or pressures). This can be enabled or disabled using the eleventh parameter in the simulation setup.

## Visualization

After running the simulation, you can visualize the results using the provided Python script `1.py`.

### Running the Visualization Script

Execute the following command:

```bash
python3 1.py
```

### Generated Plots

The script generates the following visualizations:

1. **Animation of the Entire Simulation**: Shows the evolution of a selected variable over time. You can change the variable to be plotted by modifying line 54 in `1.py`.
2. **Final Time Slice - ρ vs. x Plot**: A line plot of the density (`rho`) against the x-position at the final time step.
3. **Final Time Slice - ρ Contour Plot**: A contour plot showing the spatial distribution of density at the final time step.

### Customizing Plots

- **Changing Variables**: To plot different variables (e.g., pressure, velocity components), modify the variable selection in the script.
- **Plot Settings**: You can adjust plot titles, labels, color maps, and other aesthetics within the script to suit your preferences.

## File Structure and Code Organization

- `conservationform.C`: Defines the conservation form of the equations.
- `initstate.C`: Contains functions for initializing the state variables.
- `Solver.C`: Base class for solvers, containing common functions and interfaces.
- `Solver_HLLC.C`: Implementation of the HLLC solver.
- `Parallel.C`: Contains functions and definitions for parallel execution using OpenMP.
- `main.C`: The main function where the simulation is configured and executed.
- `1.sh`: Shell script for building and running the simulation.
- `Makefile`: Instructions for compiling the code.
- `1.py`: Python script for visualizing the simulation results.

## Usage Example

1. **Build the Project**:

   ```bash
   ./1.sh
   ```

2. **Edit Simulation Parameters**:

   Open `main.C` or the solver classes and adjust parameters as needed.

3. **Run the Simulation**:

   ```bash
   ./bin/1
   ```

4. **Visualize Results**:

   ```bash
   python3 1.py
   ```

5. **Modify and Re-run**:

   Tweak parameters, initial conditions, boundary conditions, or limiter functions, rebuild, and rerun the simulation to observe different behaviors.
