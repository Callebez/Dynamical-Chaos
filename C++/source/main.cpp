#include<iostream>
#include<fstream>
#include<iterator> // for iterators
#include<vector> // for vectors
#include<cmath>
// #include "../include/gnuplot-iostream.h"
#include "../include/plotting.hpp"
#include "../include/biffurcation.hpp"
#include "../include/penduli.hpp"
#include "../include/plotting.hpp"

#include "../include/LyapExp.hpp"
std::vector<double> lorenz(std::vector<double> coord, double rho)
{
    std::vector<double> coord_dot(3,0);
    coord_dot[0] = 10.0*(coord[1]- coord[0]);
    coord_dot[1] = coord[0]*rho- coord[0]*coord[2] - coord[1];
    coord_dot[2] = coord[0]*coord[1] - 8.0/3.0*coord[2];
    return coord_dot;
    
}
#include "../include/LyapExp.hpp"
#include "../include/lorenz.hpp"
#include "../include/printing.hpp"



void allInOne(std::vector<double>(*function)(std::vector<double>, double),
                                        std::vector<std::vector<double>> (*jacobian)(std::vector<double>&,double),
                                        double paramRange[2], std::vector<double> initialCond, 
                                        double paramStep, double integrationStep, int systemDimension)
{

    int paramIterations = (int)((fabs(paramRange[1]-paramRange[0])/paramStep));
    int integrationIterations = (int)((double)100.0/integrationStep);
    std::vector<std::vector<double>> biffurcation (4);
    std::vector<std::vector<double>> auxCoord(integrationIterations);
    auxCoord[0]= initialCond;

    std::vector<double> param (paramIterations,0);
    std::vector<double>::iterator paramValue;
    std::vector<double> xCoord (integrationIterations-1,0);
    std::vector<double> auxIntegration;

    std::vector<std::vector<double>> w = identityMatrix(initialCond.size());  
    std::vector<std::vector<double>> J = identityMatrix(initialCond.size());
    std::vector<std::vector<long double>> lyapunovExponents (paramIterations ,std::vector<long double> (initialCond.size()+1,0));
    
    for(int i = 0; i < paramIterations; i++)
    {
        param[i] = i*paramStep;
        // std::cout<<param[i]<<" \n";
    } 
    
    int k = 0;
    for(paramValue = param.begin(); paramValue<param.end(); paramValue++)
    {
        //integration for transient state
        for(int i = 0; i < integrationIterations-1; i++)
        {
            auxCoord[i+1] = rungeKutta4thSquare(function, auxCoord[i], *paramValue, integrationStep,systemDimension);
            J = jacobian(auxCoord[i],integrationStep);
            w = matMult(J,w);
            transpostSquare(w);
            gramSchmidt(w);
        
            for(uint l = 0; l < initialCond.size(); l++)
            {
                lyapunovExponents[k][l+1] +=  log(normOf(w[l]));
            }
        
            gramSchmidtNormal(w);
            transpostSquare(w);
            xCoord[i] = auxCoord[i][2];
        }
        lyapunovExponents[k][0] =(double)(*paramValue);
        for(uint j = 1; j < initialCond.size(); j++)
        {
            lyapunovExponents[k][j] = lyapunovExponents[k][j]/((long double)(100.0));
        }
        k++;
        auxCoord[0] = initialCond;
        // auxCoord[0] = auxCoord[integrationIterations-1];
        

        // for(int i = 1; i < integrationIterations-2; i++)
        // {
        //     if((xCoord[i-1]<xCoord[i]) & (xCoord[i]>xCoord[i+1]))
        //     {
        //         biffurcation[1].push_back(xCoord[i]);
        //         biffurcation[0].push_back((double)(*paramValue));
        //         if(xCoord[i+1]<xCoord[i+2])
        //         {
        //             biffurcation[2].push_back(xCoord[i+1]);
        //             biffurcation[3].push_back((double)(*(paramValue + 1)));
        //             i++;
        //         }
        //     }
        //     else if((xCoord[i]<xCoord[i-1]) & (xCoord[i]< xCoord[i+1]))
        //     {
        //         biffurcation[3].push_back(xCoord[i]);
        //         biffurcation[2].push_back((double)(*paramValue));
        //         if(xCoord[i+1]>xCoord[i+2])
        //         {
        //             biffurcation[0].push_back(xCoord[i+1]);
        //             biffurcation[1].push_back((double)(*(paramValue + 1)));
        //             i++;
        //         }
        //     }
        // }
        xCoord.clear();
    }
  

    printLyapunovsToFile(lyapunovExponents, "lyapunovExpHistLorenz1");
    printBiffucationToFile(biffurcation, "biffurcationAndLyapunovExponents1");
    plotBiffucationAndLyapunovExp("biffurcationAndLyapunovExponents1", "lyapunovExpHistLorenz1",\
                                 "Lorenz System", "{/Symbol r}", "");


}
int main()
{
    std::vector<double> integrationAux = {10.0,1.0,0.0};
    double time[2] = {0,10.0};
    // std::vector<long double>  lya = lyapunovSpectrum(lorenz,lorenzJacobian,integrationAux,0.001, 28.0);
    // std::cout<<"lyapunov numbers"<<lya[0]<<", "<<lya[1]<<", "<<lya[2]<<"\n";
    // std::cout<<"lyapunov exponents"<<exp(lya[0])<<", "<<exp(lya[1])<<", "<<exp(lya[2]);

    allInOne(lorenz,lorenzJacobian,time, integrationAux,1,0.01,3);
    // std::vector<std::vector<double>> id = biffurcation(lorenz,time,integrationAux,0.1,0.01,3);

    // printBiffucationToFile(id,"biffucationLorenz");
    // plotBiffucation("biffucationLorenz","Lorenz System", " {/Symbol r}", " ");
    // plot2D("./outputs/images/biffucationLorenz","./outputs/txt/biffucationLorenzmax.dat", "max points"," plot ./outputs/txt/biffucationLorenzmin.dat w dots title \'min points\'  set title \'Biffurcation Diagram for the Lorenz system\'" );
    
    return 0;    
}
