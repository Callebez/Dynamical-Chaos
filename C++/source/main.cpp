#include<iostream>
#include<fstream>
#include<iterator> // for iterators
#include<vector> // for vectors
#include<cmath>
#include<string>
// #include "../include/gnuplot-iostream.h"
 //#include "../include/rungekutta4thSquare.hpp"
#include "../include/biffurcation.hpp"
#include "../include/penduli.hpp"
#include "../include/LyapExp.hpp"
#include "../include/lorenz.hpp"
#include "../include/LinearAlgebra.hpp"

int main()
{

    // std::vector<std::vector<double>> M ={{2, 4}, 
                                        // {3, 4}};
    std::vector<double> integrationAux = {1.0,2.0,3.0};

    std::vector<long double> lya = lyapunovSpectrum(lorenz,lorenzJacobian,integrationAux,0.001,28.0);
    // std::vector<std::vector<double>> A = matMult(M,N);
    std::cout<<exp(lya[0])<<", "<<exp(lya[1])<<", "<<exp(lya[2])<<"\n ";
    // std::cout<<"sum:"<< lya[0]
    // std::vector<std::vector<double>> m = identityMatrix(5);
    // for(int i = 0; i < 3; i++)
    // {s
    
    // }
    // printMatrix(N);
  
    // std::vector<double> integrationAux = {1,1,1};
    // std::vector<double> integrationAux2 = {20,4,0};
    // normalize(integrationAux2);
    // std::cout<< integrationAux2[0]<<", "<< integrationAux2[1] << ", "<< integrationAux2[2];

    // double time[2] = {0,10000.0};
    // completeRungeKuttaToFile(classicalPendulum,integrationAux, 1,0.01,4,time);
    // std::cout<<lyapunov(lorenz,integrationAux,lorenzJaconian,3,0.01);
    /*
    // std::vector<double> auxVec (4,0);
    // double step = 1e-3;
    
    // double timespan[] = {0,10};
    // auxVec = rungeKutta4thSquare(classicalPendulum, integrationAux, 1e-3, step, 4);

    // completeRungeKuttaToFile(classicalPendulum, integrationAux, 1, step, 4, timespan);


    // std::vector<double> param (10,0);
    // std::vector<double>::iterator paramValue;
    // for(int i = 0; i < 10; i++)
    // {
    //     param[i] = i+0.1;
    // } 
    // for(paramValue = param.begin(); paramValue<param.end(); paramValue++)
    // {
    //     std::cout<<"atual: "<<*paramValue<<"\n";
    //     std::cout<<"seguinte: "<<*(paramValue+1)<<"\n";
    // }

//     std::vector<std::vector<double>> Matrix;  
//     std::vector<double> line (4, 0);
//     std::vector<double> line2 (4, 1);
//     std::vector<double> aux;
//     double x;
//     Matrix.push_back(line);
//     Matrix.push_back(line2);
//     Matrix[0].push_back(1);
//     Matrix[1].push_back(0);
//     for(uint i = 0; i < Matrix.size(); i++)
//     {
//         for(uint j = 0; j < Matrix[0].size(); j++)
//         {
//             std::cout<<Matrix[i][j] <<" ";
//         }
//         std::cout<<std::endl;
//     }
    
//     aux = Matrix[1];
//     Matrix.clear();
//     Matrix.push_back(aux);
// std::cout<<std::endl;
//     for(uint i = 0; i < Matrix.size(); i++)
//     {
//         for(uint j = 0; j < Matrix[0].size(); j++)
//         {
//             std::cout<<Matrix[i][j] <<" ";
//         }
//         std::cout<<std::endl;
//     }
    double range[2] = {0.0 , 2.0};
    std::vector<std::vector<double>> biff = biffurcation(classicalPendulum,range,integrationAux, 0.001,0.1,4);
    // int iterations = (int)(abs(time_span[1]-time_span[0])/step);
    std::ofstream biffdiagramA;
    std::ofstream biffdiagramB;
    biffdiagramA.open("DAT/biffurcationsC.dat");
    biffdiagramB.open("DAT/biffurcationsD.dat");
    //std::vector<double> auxVec = rungeKutta4thSquare(function, initialCond, param, step, dimension);
    for(int j = 0; j < (int)biff[0].size(); j++)
    {
        param[i] = i*0.1;
        // std::cout<<param[i]<<" \n";
    } 
    double time_span[2] = {0.0,30.0};
    double step = 0.0001;
    int iterations = (int)(fabs(time_span[1]-time_span[0])/step);
    std::vector<double>time(iterations);
    for(int i = 0; i < iterations; i++)
    {   
        time[i] = i*step;

    }
    std::ofstream Xt;
    std::ofstream Yt;
    std::ofstream Pxt;
    std::ofstream Pyt;
    Xt.open("classicalXt.dat");
    Yt.open("classicalYt.dat");
    Pxt.open("classicalPxt.dat");
    Pyt.open("classicalPyt.dat");

     for(paramValue = param.begin(); paramValue<param.end(); paramValue++)
    {
        integrationAux = {6,0,1,1};
        for(int j = 0; j < iterations; j++)
        {
            integrationAux = rungeKutta4thSquare(classicalPendulum, integrationAux, *paramValue, step, 4);
            // for(int i = 0; i < 4; i++)
            // {
            Xt << time[j]<<"  ";
            Xt << integrationAux[0] <<"  ";
            
            Pxt << time[j]<<"  ";
            Pxt << integrationAux[1] <<"  ";

            Yt << time[j]<<"  ";
            Yt << integrationAux[2] <<"  ";
            
            Pyt << time[j]<<"  ";
            Pyt << integrationAux[3] <<"  ";
            // }
            Xt<< std::endl;
            Yt<< std::endl;
            Pxt<< std::endl;
            Pyt<< std::endl;
        }
        Xt<<"\n" <<std::endl;
        Yt<<"\n" <<std::endl;
        Pxt<<"\n" <<std::endl;
        Pyt<<"\n" <<std::endl;
    }
    biffdiagramA.close();
    biffdiagramB.close();*/
    // int iterations=1e4;
    // double step = 1e-3;
    // std::vector<double> coord=rungeKutta4thSquare(classicalPendulum, integrationAux, 1e-3, step, 4);
    // std::vector<double> Exponents=LyapunovExponents(classicalPendulum,integrationAux,step,iterations);
    
    return 0;    
}
