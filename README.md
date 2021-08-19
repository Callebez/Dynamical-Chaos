# Dynamical Chaos
## A project for dynamical exploration of chaos

Dynamical Chos is a project developed with the support of CNPq, made as a one-year project for PIBIC-UnB (research internship at the University of Brasília) carried out by Callebe Reis and Shohei Sato under the guidance of Professor Olavo Leopoldino da Silva Filho.

The main goal of the project is to provide a set of tools to analyze autonomous dynamic systems and the presence of chaos. Having tools for numeric integration of ordinary differential equations (RK45) and (RK4), finding the Lyapunov exponent spectrum, creating bifurcation diagrams and plotting 2D and 3D functions using a gnuplot pipe, the project achieves its goals, being able to provide valuable information about the behavior of various dynamic systems.

## Getting started 
The only dependencies are [gnuplot](http://www.gnuplot.info/) and [make](https://www.gnu.org/software/make/). 

<ins>**1. Downloading the repository:**</ins>

Clone the repository in the desired location with `git clone https://github.com/CallebeRReis/DynamicalSystems.git`.

<ins>**2. Defining the dynamical system of interest:**</ins>

It is possible that the system you're looking for isn't implemented yet, in that case you must add the system by creating a `source/systemname.cpp` and a  `include/systemname.hpp files`. 
The structure of the function must be something like: 

```
std::vector<double> systemName(std::vector<double> coordinates, double parameter)
{
    std::vector<double> dot_coordinates (size of your system)

    dot_coordinates[0] = ... 
    dot_coordinates[1] = ...
    dot_coordinates[2] = ...
    .
    .
    . 
    dot_coordinates[size of your system] = ...
    return dot_coordinates
}
```
if you want to calculate the Lyapunov expoents spectrum, you may want to include a function with the jacobian * step + identity matrix, where the identity matrix have the same size as the jacobian:

```
std::vector<std::vector<double>> systemJacobian(std::vector<double>& coordinates, double parameter, double step)
{
    std::vector<std::vector<double>> jacobian (size of your system, std::vector<double>(size of your system,0));

    jacobian[0][0] =  equation*step + 1; 
    jacobian[0][1] =  ...
    .
    .
    .
    jacobian[1][0] = ...
    jacobian[1][1] = equation * step + 1;
    jacobian[1][2] = ...
    .
    .
    .
    jacobian[size of your system][size of your system] = equation * step + 1;
    .
    .
    .
    
    return jacobian;
}
```
<ins>**3. Making the analysis:**</ins>

In the `source/main.cpp` include the `#include "../include/systemname.hpp"` declare the initial conditions by `std::vector<double> initialConditions = {a,b,...,c}` where "a", "b" and "c" are numbers.
 
declare the space in wich the parameter will vary `double  Paramter Interval [2] = {a, b}`, 
 
<ins>**3.1 Simple plot:**</ins>
declare a matrix for storing the integration data, and a vector for the time variable
```
  std::vector<std::vector<double>> integration Data
  std::vector<double> Time
```
then 
```
  completeRungeKutta45(Function in wich the system is describe,
                       Initial Conditions, 
                       Integration Data, 
                       Parameter Value, 
                       Initial Step, 
                       System Dimension,
                       Integration Tolerance, 
                       TimeInterval,
                       Iime);
                       
  printMatrixToFile(Integration Data,
                    "name_of_the_file.dat");
                    
  plot3D("output_File_Name",
        "name_of_the_file.dat", 
        "Name of the System", 
        "optinal gnuplot parameters");
```
Then, open the terminal in the folder of the program, and type `make run`, the image should be in `outputs/output_File_Name.png`
 
<ins>**3.2 bifurcation plot:**</ins>
The bifurcation plot will depend on the parameter interval. 
``` 
    std::vector<std::vector<double>> Bifurcation Data = bifurcation(Function in wich the system is describe,
                                                                    Paramter Interval,
                                                                    Initial Conditions,
                                                                    Step of the variation of the paramter,
                                                                    Coordinate beeing analysed,
                                                                    System dimension);

    printBiffucationToFile(Bifurcation Data,
                          "systems_bifucation_file_name");
    plotBiffucation("systems_bifucation_file_name",
                    "Name of the System",
                    " Name of the parameter varying",
                    "Gnuplot optinal paramters");
``` 
<ins>**3.2 Lyapunov Exponent Spectrum:**</ins>
for calculating the spectrum, one must do:
```
    std::vector<long double> Lyapunov Spectrum = lyapunovSpectrum(Function in wich the system is describe,
                                                                  Jacobian of the system,
                                                                  Initial Conditions,
                                                                  Integration tolerance,
                                                                  Time elapsed,
                                                                  Initial Step,
                                                                  Function parameter);      
```
and for printing the results:
```
    std::cout << " Lyapunov numbers: " << Lyapunov Spectrum[0] << ", " << Lyapunov Spectrum[1] << ", " << Lyapunov Spectrum[2]<< "\n";
    std::cout << " Lyapunov exponents: " << exp(Lyapunov Spectrum[0]) << ", " << exp(Lyapunov Spectrum[1]) << ", " << exp(Lyapunov Spectrum[2]) << "\n";

``

Feel free to send us aditional systems that aren't implemented yet!
