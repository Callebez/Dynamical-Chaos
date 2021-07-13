#include<iostream>
#include<fstream>
#include<iterator> // for iterators
#include<vector> // for vectors
#include<cmath>
// #include "../include/gnuplot-iostream.h"
// #include "../include/rungekutta4thSquare.hpp"
#include "../include/biffurcation.hpp"
#include "../include/penduli.hpp"
#include "../include/LyapExp.hpp"
int main()
{

    std::vector<double> integrationAux = {6,1,0,1};
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
      
        for(int i = 0; i < 2; i++)
        {
            biffdiagramA << biff[i][j] <<"  ";
            biffdiagramB << biff[2+i][j] <<"  ";
        }
        biffdiagramA<< std::endl;
        biffdiagramB<< std::endl;
    }
    biffdiagramA.close();
    biffdiagramB.close();*/
    int iterations=1e4;
    double step = 1e-3;
    std::vector<double> coord=rungeKutta4thSquare(classicalPendulum, integrationAux, 1e-3, step, 4);
    std::vector<double> Exponents=LyapunovExponents(classicalPendulum,integrationAux,step,iterations);
    
    return 0;    
}
