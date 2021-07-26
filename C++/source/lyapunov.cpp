
#include "../include/LyapExp.hpp"
#include "../include/LinearAlgebra.hpp"
/* 
    lyapunovSpectrum:   Function for the obtain the lyapunov spectrum of exponents
    function:           Set of differential equations
    jacobian:           Jacobian matrix of the system
    initialCond:        Initial conditions of the system. Similiar initial conditions may have different spectruns
    step:               Step size of the integration
    param:              Paramater of interrest when seeing the how the system biffurcates

 
*/
std::vector<long double> lyapunovSpectrum(std::vector<double> (*function)(std::vector<double>, double),
                                std::vector<std::vector<double>> (*jacobian)(std::vector<double>&,double), 
                                std::vector<double>& initialCond, double step, double param)
{
   
    uint iterations = (uint)(100.0/step);
    std::vector<double> coord = initialCond;
    std::vector<long double> lyapunovExponents (initialCond.size(),0);

    std::vector<std::vector<double>> w = identityMatrix(initialCond.size());  
    std::vector<std::vector<double>> J = identityMatrix(initialCond.size());
    
    //loop for making the initial conditions sit near the atractor  
    for(uint i = 0; i < 1000; i++)
    {
        coord = rungeKutta4thSquare(function, coord, param, 0.01, initialCond.size());
    }
    // main loop
    for(uint i = 0; i < iterations; i++)
    {
        //runge kutta does the evolution of the system 
        coord = rungeKutta4thSquare(function, coord, param, step, initialCond.size());
        //the jacobian times the matrix makes the vectors in w align with the directions 
        //of the eigenvectors  
        J = jacobian(coord,step);
        w = matMult(J,w);

        transpostSquare(w);
        gramSchmidt(w);

        for(uint j = 0; j < initialCond.size(); j++)
        {
            lyapunovExponents[j] +=  log(normOf(w[j]));
        }
    
        gramSchmidtNormal(w);
        transpostSquare(w);
    }
    for(uint j = 0; j < initialCond.size(); j++)
    {
        lyapunovExponents[j] = lyapunovExponents[j]/((long double)(100.0));
    }
    
    return lyapunovExponents;
}