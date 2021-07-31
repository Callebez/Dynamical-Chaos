
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
                                std::vector<std::vector<double>> (*jacobian)(std::vector<double>&, double,double), 
                                std::vector<double>& initialCond, double tol, double time, 
                                double step, double param)
{
   
    std::vector<double> coord = initialCond;
    std::vector<long double> lyapunovExponents (initialCond.size(),0);
    std::vector<double> hs;

    std::vector<std::vector<double>> w = identityMatrix(initialCond.size());  
    std::vector<std::vector<double>> J = identityMatrix(initialCond.size());

    std::vector<std::vector<double>> rk45;
    completeRungeKutta45(function,initialCond, rk45, param,step,coord.size(),tol, time, hs);
    initialCond = rk45.back();
    // main loop
    for(uint i = 0; i < rk45.size(); i++)
    {
        
        //runge kutta does the evolution of the system 
        coord = rk45[i];
        //the jacobian times the matrix makes the vectors in w align with the directions 
        //of the eigenvectors  
        J = jacobian(coord, param,hs[i]);
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
void laypunovVaringParameter(std::vector<double>(*function)(std::vector<double>, double),
                                            std::vector<std::vector<double>> (*jacobian)(std::vector<double>&, double,double), 
                                            double paramRange[2], std::vector<double> initialCond, double tol,
                                            double paramStep, int coordBeingAnalysed, int systemDimension,
                                            std::vector<std::vector<double>>& lyapunovSpectrumRange)
{

    int paramIterations = (int)((fabs(paramRange[1]-paramRange[0])/paramStep));    

   //inicialization for the integration

    std::vector<double> param (paramIterations,0);
    std::vector<double>::iterator paramValue;
    std::vector<double> xCoord;
    double integrationStep = 0.1;
    for(int i = 0; i < paramIterations; i++)
    {
        param[i] = paramRange[0] + i*paramStep;
    } 
    std::vector<std::vector<double>> rk45;
    std::vector<double> hs;
    std::vector<std::vector<long double>> lyapunovExponents ;

    lyapunovSpectrumRange.resize(param.size(),std::vector<double>(initialCond.size()+1));
    for(paramValue = param.begin(); paramValue<param.end(); paramValue++)
    {
        lyapunovExponents.emplace_back( lyapunovSpectrum(function,jacobian,initialCond,tol,10,integrationStep,*paramValue));
    }
    for(int i = 0; i < param.size();i++)
    {
        lyapunovSpectrumRange[i][0] = param[i];

        for(int j = 1; j < (int)initialCond.size() +1; j++)
        {
            lyapunovSpectrumRange[i][j] = lyapunovExponents[i][j-1];
        }
    }
}