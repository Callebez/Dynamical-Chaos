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



// void allInOne(std::vector<double>(*function)(std::vector<double>, double),
//                                         std::vector<std::vector<double>> (*jacobian)(std::vector<double>&,double),
//                                         double paramRange[2], std::vector<double> initialCond, 
//                                         double paramStep, double integrationStep, int systemDimension)
// {

//     int paramIterations = (int)((fabs(paramRange[1]-paramRange[0])/paramStep));
//     int integrationIterations = (int)((double)100.0/integrationStep);
//     std::vector<std::vector<double>> biffurcation (4);
//     std::vector<std::vector<double>> auxCoord(integrationIterations);
//     auxCoord[0]= initialCond;

//     std::vector<double> param (paramIterations,0);
//     std::vector<double>::iterator paramValue;
//     std::vector<double> xCoord (integrationIterations-1,0);
//     std::vector<double> auxIntegration;

//     std::vector<std::vector<double>> w = identityMatrix(initialCond.size());  
//     std::vector<std::vector<double>> J = identityMatrix(initialCond.size());
//     std::vector<std::vector<long double>> lyapunovExponents (paramIterations ,std::vector<long double> (initialCond.size()+1,0));
    
//     for(int i = 0; i < paramIterations; i++)
//     {
//         param[i] = i*paramStep;
//         // std::cout<<param[i]<<" \n";
//     } 
    
//     int k = 0;
//     for(paramValue = param.begin(); paramValue<param.end(); paramValue++)
//     {
//         //integration for transient state
//         for(int i = 0; i < integrationIterations-1; i++)
//         {
//             auxCoord[i+1] = rungeKutta4thSquare(function, auxCoord[i], *paramValue, integrationStep,systemDimension);
//             J = jacobian(auxCoord[i],integrationStep);
//             w = matMult(J,w);
//             transpostSquare(w);
//             gramSchmidt(w);
        
//             for(uint l = 0; l < initialCond.size(); l++)
//             {
//                 lyapunovExponents[k][l+1] +=  log(normOf(w[l]));
//             }
        
//             gramSchmidtNormal(w);
//             transpostSquare(w);
//             xCoord[i] = auxCoord[i][2];
//         }
//         lyapunovExponents[k][0] =(double)(*paramValue);
//         for(uint j = 1; j < initialCond.size(); j++)
//         {
//             lyapunovExponents[k][j] = lyapunovExponents[k][j]/((long double)(100.0));
//         }
//         k++;
//         auxCoord[0] = initialCond;
//         // auxCoord[0] = auxCoord[integrationIterations-1];
        

//         // for(int i = 1; i < integrationIterations-2; i++)
//         // {
//         //     if((xCoord[i-1]<xCoord[i]) & (xCoord[i]>xCoord[i+1]))
//         //     {
//         //         biffurcation[1].push_back(xCoord[i]);
//         //         biffurcation[0].push_back((double)(*paramValue));
//         //         if(xCoord[i+1]<xCoord[i+2])
//         //         {
//         //             biffurcation[2].push_back(xCoord[i+1]);
//         //             biffurcation[3].push_back((double)(*(paramValue + 1)));
//         //             i++;
//         //         }
//         //     }
//         //     else if((xCoord[i]<xCoord[i-1]) & (xCoord[i]< xCoord[i+1]))
//         //     {
//         //         biffurcation[3].push_back(xCoord[i]);
//         //         biffurcation[2].push_back((double)(*paramValue));
//         //         if(xCoord[i+1]>xCoord[i+2])
//         //         {
//         //             biffurcation[0].push_back(xCoord[i+1]);
//         //             biffurcation[1].push_back((double)(*(paramValue + 1)));
//         //             i++;
//         //         }
//         //     }
//         // }
//         xCoord.clear();
//     }
  

//     printLyapunovsToFile(lyapunovExponents, "lyapunovExpHistLorenz1");
//     printBiffucationToFile(biffurcation, "biffurcationAndLyapunovExponents1");
//     plotBiffucationAndLyapunovExp("biffurcationAndLyapunovExponents1", "lyapunovExpHistLorenz1",
//                                  "Lorenz System", "{/Symbol r}", "");


// }
int main()
{
    std::vector<double> integrationAux = {1.0,2.0,3.0};
    // double time[2] = {0,10.0};
    // std::vector<std::vector<double>> rk45 (100000, std::vector<double>(3,0));
    // double time;
    std::vector<double> hs; 
    // std::vector<std::vector<double>>  biff = biffurcation(lorenz, time,integrationAux,0.1,0.1,0,3);
    // printBiffucationToFile(biff, "testplot");
    // plotBiffucation("testplot","lorenz test","rho", "");
    // completeRungeKutta45(lorenz,integrationAux,100000, rk45, 28.0,0.001,3,1e-5, time, hs);
    // printMatrixToFile(rk45,"testeRK45.dat");
    // plot3D("testeRK45","testeRK45", "lorenz Atrractor", "");

  //  printMatrix(rk45);
    // double stepNew; 
    // rungeKutta45(lorenz,integrationAux,rk45,28.0, 0.0150955,3, 1e-2 ,stepNew);
    // std::vector<double> rk4s =rungeKutta4thSquare(lorenz,integrationAux,28.0, 1e-2,3);
    // std::cout<<"rk4 novo: "<< rk45[0]<<", "<<rk45[1]<<", "<<rk45[2]<<"\n";
    // std::cout<<"error : "<< rk5[0]<<", "<<rk5[1]<<", "<<rk5[2]<<"\n";
    // std::cout<<"rk4 velho: "<< rk4s[0]<<", "<<rk4s[1]<<", "<<rk4s[2]<<"\n";
    // std::cout<<"novo step: "<<stepNew;
    double time[2] = {0.0,100.0};
    // int iterations = (int)(fabs(time[1]-time[0])/0.1);
    std::vector<std::vector<double>> A;// (iterations, std::vector<double> (integrationAux.size()+1));
    laypunovVaringParameter(lorenz,lorenzJacobian,time,integrationAux,1e-6,0.5,0,3,A);
    printMatrixToFile(A,"arquivoTesteLyapunovVsRho.dat");

    // std::vector<long double>  lya = lyapunovSpectrum(lorenz,lorenzJacobian,integrationAux,1e-6,100,0.01, 28.0);
    // std::cout<<"lyapunov numbers: "<<lya[0]<<", "<<lya[1]<<", "<<lya[2]<<"\n";
    // std::cout<<"lyapunov exponents: "<<exp(lya[0])<<", "<<exp(lya[1])<<", "<<exp(lya[2])<<"\n";
    // std::cout<<"their sum : "<< lya[0]+lya[1]+lya[2];
  
    // std::vector<std::vector<double>> biff = biffurcation(lorenz, time, integrationAux, 0.1,0.001,2,3);

    // printBiffucationToFile(biff,"biffucationLorenzTesteRK45");
    // plotBiffucation("biffucationLorenzTesteRK45","Lorenz System", " {/Symbol r}", " ");
    // plot2D("./outputs/images/biffucationLorenz","./outputs/txt/biffucationLorenzmax.dat", "max points"," plot ./outputs/txt/biffucationLorenzmin.dat w dots title \'min points\'  set title \'Biffurcation Diagram for the Lorenz system\'" );
    
    return 0;    
}
