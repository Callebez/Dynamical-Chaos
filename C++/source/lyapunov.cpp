
#include "../include/LyapExp.hpp"
#include "../include/LinearAlgebra.hpp"
#include "../include/printing.hpp"
#include "../include/plotting.hpp"

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
                                std::vector<double>& initialCond, double tol, double time, 
                                double step, double param)
{
   
    std::vector<double> coord = initialCond;
    std::vector<long double> lyapunovExponents (initialCond.size(),0);
    std::vector<double> hs;

    std::vector<std::vector<double>> w = identityMatrix(initialCond.size());  
    std::vector<std::vector<double>> J = identityMatrix(initialCond.size());

    std::vector<std::vector<double>> rk45;
    completeRungeKutta45(function,initialCond,0, rk45, param,step,coord.size(),tol, time, hs);
    
    // main loop
    for(uint i = 0; i < rk45.size(); i++)
    {
        
        //runge kutta does the evolution of the system 
        coord = rk45[i];
        //the jacobian times the matrix makes the vectors in w align with the directions 
        //of the eigenvectors  
        J = jacobian(coord,hs[i]);
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
        lyapunovExponents[j] = lyapunovExponents[j]/(time);
    }
    return lyapunovExponents;
}