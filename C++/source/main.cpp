#include<iostream>
#include<fstream>
#include<iterator> // for iterators
#include<vector> // for vectors
#include<cmath>
// #include "../include/gnuplot-iostream.h"
// #include "../include/rungekutta4thSquare.hpp"
#include "../include/biffurcation.hpp"
#include "../include/penduli.hpp"
std::vector<double> lorenz(std::vector<double> coord, double rho)
{
    std::vector<double> coord_dot(3,0);
    coord_dot[0] = 10.0*(coord[1]- coord[0]);
    coord_dot[1] = coord[0]*rho- coord[0]*coord[2] - coord[1];
    coord_dot[2] = coord[0]*coord[1] - 8.0/3.0*coord[2];
    return coord_dot;
    
}
int main()
{

    std::vector<double> integrationAux = {1,1,1};
    // std::vector<double> auxVec (4,0);
    // double step = 1e-3;
    
    // double timespan[] = {0,10};
    // auxVec = rungeKutta4thSquare(classicalPendulum, integrationAux, 1e-3, step, 4);
    // completeRungeKuttaToFile(lorenz, integrationAux, 1, step, 4, timespan);


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
    
    std::vector<double> param (400,0);
    std::vector<double>::iterator paramValue;

    for(int i = 0; i < 400; i++)
    {
        param[i] = i*0.5;
        // std::cout<<param[i]<<" \n";
    } 
    double time_span[2] = {0.0,50.0};
    int iterations = (int)(fabs(time_span[1]-time_span[0])/0.001);
    std::ofstream fileName;
    fileName.open("outputrk4thnew.dat");

     for(paramValue = param.begin(); paramValue<param.end(); paramValue++)
    {
        integrationAux = {1,1,1};
        for(int j = 0; j < iterations; j++)
        {
            integrationAux = rungeKutta4thSquare(lorenz, integrationAux, *paramValue, 0.01, 3);
            for(int i = 0; i < 3; i++)
            {
                fileName << integrationAux[i] <<"  ";
            }
            fileName<< std::endl;
        }
        fileName<<"\n" <<std::endl;
    }
  

    fileName.close();

   //completeRungeKuttaToFile(lorenz, {1,1,1}, 28, 0.01, 3, time_span);
    // double range[2] = {0 , 20};
    // std::vector<std::vector<double>> biff = biffurcation(classicalPendulum,range,integrationAux, 0.01,0.01,4);
    // // int iterations = (int)(abs(time_span[1]-time_span[0])/step);
    // std::ofstream biffdiagramA;
    // std::ofstream biffdiagramB;
    // biffdiagramA.open("biffurcationsC.dat");
    // biffdiagramB.open("biffurcationsD.dat");
    // //std::vector<double> auxVec = rungeKutta4thSquare(function, initialCond, param, step, dimension);
    // for(int j = 0; j < (int)biff[0].size(); j++)
    // {
      
    //     for(int i = 0; i < 2; i++)
    //     {
    //         biffdiagramA << biff[i][j] <<"  ";
    //         biffdiagramB << biff[2+i][j] <<"  ";
    //     }
    //     biffdiagramA<< std::endl;
    //     biffdiagramB<< std::endl;
    // }
    // biffdiagramA.close();
    // biffdiagramB.close();
    return 0;    
}
