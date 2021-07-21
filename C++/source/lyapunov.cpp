
#include "../include/LyapExp.hpp"
#include "../include/LinearAlgebra.hpp"

std::vector<long double> lyapunovSpectrum(std::vector<double> (*function)(std::vector<double>, double),
                                std::vector<std::vector<double>> (*jacobian)(std::vector<double>&,double), 
                                std::vector<double>& initalCond, double step, double param)
{
    // double step = 0.001;
    uint iterations = (uint)(100.0/step);
    std::vector<double> coord = initalCond;
    std::vector<long double> lyapunovExponents (initalCond.size(),0);

    std::vector<std::vector<double>> w = identityMatrix(initalCond.size());  
    std::vector<std::vector<double>> J = identityMatrix(initalCond.size());
    for(uint i = 0; i < 100; i++)
    {
        coord = rungeKutta4thSquare(function, coord, param, step, initalCond.size());
    }

    for(uint i = 0; i < iterations; i++)
    {
        coord = rungeKutta4thSquare(function, coord, param, step, initalCond.size());
        J = jacobian(coord,step);
        w = matMult(J,w);
        transpostSquare(w);
        gramSchmidt(w);
    
        for(uint j = 0; j < initalCond.size(); j++)
        {
            lyapunovExponents[j] +=  log(normOf(w[j]));
        }
    
        gramSchmidtNormal(w);
        transpostSquare(w);
    }
    for(uint j = 0; j < initalCond.size(); j++)
    {
        lyapunovExponents[j] = lyapunovExponents[j]/((long double)(100.0));
    }
    
    return lyapunovExponents;
}