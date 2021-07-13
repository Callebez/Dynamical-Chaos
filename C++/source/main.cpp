#include<iostream>
#include<fstream>
#include<iterator> // for iterators
#include<vector> // for vectors
#include<cmath>
#include<string>
// #include "../include/gnuplot-iostream.h"
// #include "../include/rungekutta4thSquare.hpp"
#include "../include/biffurcation.hpp"
#include "../include/penduli.hpp"
#include "../include/plotting.hpp"

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
    
    std::vector<double> param (20,0);
    std::vector<double>::iterator paramValue;

    for(int i = 0; i < 20; i++)
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
    plotAnimate2D("classicalXt.gif", "classicalXt.dat", 30);
    plotAnimate2D("classicalYt.gif", "classicalYt.dat", 30);
    plotAnimate2D("classicalPxt.gif", "classicalPxt.dat", 30);
    plotAnimate2D("classicalPyt.gif", "classicalPyt.dat", 30);
  


    Xt.close();
    Yt.close();
    Pxt.close();
    Pyt.close();

    // plotAnimate("lorenz.gif", "outputrk4thnew.dat", 10);
    // std::cout<<"success!";
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
