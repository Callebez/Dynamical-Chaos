#include "../include/rungekutta4thSquare.hpp"
#include "../include/LinearAlgebra.hpp"
#include <iostream>

std::vector<double> rungeKutta4thSquare(std::vector<double> (*function)(std::vector<double>, double), 
                                  std::vector<double> &coord, double param, double step, 
                                  int dimension )
{
    std::vector<double> k1 (dimension);
    std::vector<double> k2 (dimension);
    std::vector<double> k3 (dimension);
    std::vector<double> k4 (dimension);
    
    k1 = function(coord, param);
    std::cout<<"k1 velho: "<< k1[0]<<", "<<k1[1]<<", "<<k1[2]<<"\n";

    k2 = function(updateCoord(coord, k1, (double)step, 2.0,dimension), param);
    std::cout<<"k2 velho: "<< k2[0]<<", "<<k2[1]<<", "<<k2[2]<<"\n";

    k3 = function(updateCoord(coord, k2, (double)step, 2.0,dimension), param);
    std::cout<<"k3 velho: "<< k3[0]<<", "<<k3[1]<<", "<<k3[2]<<"\n";

    k4 = function(updateCoord(coord, k3, (double)step, 1.0,dimension), param);
    std::cout<<"k4 velho: "<< k4[0]<<", "<<k4[1]<<", "<<k4[2]<<"\n\n";

    for(int i = 0; i < dimension; i++)
    {
        coord[i] = coord[i] + ((double)step/(double)6.0)*(k1[i]+2.0*k2[i] + 2.0*k3[i]+ k4[i]);
    }
    return coord;
}
std::vector<double> rungeKutta4thSquarePertubation(std::vector<double> (*function)(std::vector<double>,std::vector<double>, double), 
                                  std::vector<double> coord,std::vector<double> pertubation, double param, double step, 
                                  int dimension, std::vector<double>& functionAval )
{
    std::vector<double> k1 (dimension);
    std::vector<double> k2 (dimension);
    std::vector<double> k3 (dimension);
    std::vector<double> k4 (dimension);
    
    k1 = function(coord, pertubation, param);

    k2 = function(updateCoord(coord, k1, (double)step, 2.0,dimension), pertubation, param);
    k3 = function(updateCoord(coord, k2, (double)step, 2.0,dimension), pertubation, param);
    k4 = function(updateCoord(coord, k3, (double)step, 1.0,dimension), pertubation, param);
    for(int i = 0; i < dimension; i++)
    {
        coord[i] = coord[i] + ((double)step/(double)6.0)*(k1[i]+2.0*k2[i] + 2.0*k3[i]+ k4[i]);
    }
    functionAval = k4;
    return coord;
}
std::vector<double> updateCoord(std::vector<double> coord, std::vector<double> increment, double step, double scale,int dimension)
{
    
    std::vector<double> updatedCoord (dimension);
    for(int i = 0; i < dimension; i++)
    {
        updatedCoord[i] = coord[i] + ((double)step/(double)scale)*(increment[i]);
    }
    return updatedCoord;
}
std::vector<double> updateCoord(std::vector<double> coord, std::vector<double> increment, double step, int dimension)
{
    
    std::vector<double> updatedCoord (dimension);
    for(int i = 0; i < dimension; i++)
    {
        updatedCoord[i] = coord[i] + step*(increment[i]);
    }
    return updatedCoord;
}
void rungeKutta45(std::vector<double> (*function)(std::vector<double>, double), 
                  std::vector<double> &coord,  double& error,
                  double param, double step, int dimension,
                  double tol, double& stepNew)
{
    // double a2 = 0.25;
    // double a3 = 0.375;
    // double a4 = 12.0/13.0;
    // double a6 = 0.5;
    double b21 = 0.25;
    double b31 = 3.0/32.0;
    double b32 = 9.0/32.0;
    double b41 = 1932.0/2197.0;
    double b42 = -7200.0/2197.0;
    double b43 = 7296.0/2197.0;
    double b51 = 439.0/216.0;
    double b52 = -8.0;
    double b53 = 3680.0/513.0;
    double b54 = -845.0/4104.0;
    double b61 = -8.0/27.0;
    double b62 = 2.0;
    double b63 = -3544.0/2565.0;
    double b64 = 1859.0/4104.0;
    double b65 = -11.0/40.0;
    double c1 = 25.0/216.0;
    double c3 = 1408.0/2565.0;
    double c4 = 2197.0/4104.0;
    double c5 = -0.20;
    double d1 = 1.0/360.0;
    double d3 = -128.0/4275.0;
    double d4 = -2197.0/75240.0;
    double d5 = 0.02;
    double d6 = 2.0/55.0;
   
    
    std::vector<double> k1 (dimension);
    std::vector<double> k2 (dimension);
    std::vector<double> k3 (dimension);
    std::vector<double> k4 (dimension);
    std::vector<double> k5 (dimension);
    std::vector<double> k6 (dimension);

    std::vector<double> aux1 (dimension);
    std::vector<double> aux2 (dimension);
    std::vector<double> aux3 (dimension);
    std::vector<double> aux4 (dimension);
    std::vector<double> aux5 (dimension);
    std::vector<double> aux6 (dimension);

    std::vector<double> coordTemp = coord;

    k1 = function(coord, param);

    k2 = function(updateCoord(coord, k1, (double)step, b21, dimension), param);

    scaleVector(k1,aux1,b31);       //aux1 = k1 * b31
    scaleVector(k2,aux2,b32);       //aux2 = k2 * b32
    
    sumVectors(aux1,aux2, aux2);    //aux2 = aux1 + aux2

    k3 = function(updateCoord(coord, aux2, (double)step,dimension), param);
    
    scaleVector(k1,aux1,b41);       //aux1 = k1 * b41
    scaleVector(k2,aux2,b42);       //aux2 = k2 * b42
    scaleVector(k3,aux3,b43);       //aux3 = k3 * b43

    sumVectors(aux1,aux2, aux2);    //aux2 = aux1 + aux2
    sumVectors(aux2,aux3, aux3);    //aux3 = aux2 + aux3 
       
    k4 = function(updateCoord(coord, aux3, (double)step, dimension), param);
    
    scaleVector(k1,aux1,b51);       //aux1 = k1 * b51
    scaleVector(k2,aux2,b52);       //aux2 = k2 * b52
    scaleVector(k3,aux3,b53);       //aux3 = k3 * b53
    scaleVector(k4,aux4,b54);       //aux4 = k4 * b54

    sumVectors(aux1,aux2, aux2);    //aux2 = aux1 + aux2
    sumVectors(aux3,aux4, aux4);    //aux4 = aux4 + aux3 
    sumVectors(aux2,aux4, aux4);    //aux4 = aux4 + aux2
    
    k5 = function(updateCoord(coord, aux4, (double)step, dimension), param);

    scaleVector(k1,aux1,b61);       //aux1 = k1 * b61
    scaleVector(k2,aux2,b62);       //aux2 = k2 * b62
    scaleVector(k3,aux3,b63);       //aux3 = k3 * b63
    scaleVector(k4,aux4,b64);       //aux4 = k3 * b64
    scaleVector(k5,aux5,b65);       //aux5 = k5 * b65

    sumVectors(aux1,aux2, aux2);    //aux2 = aux1 + aux2
    sumVectors(aux3,aux4, aux4);    //aux4 = aux4 + aux3
    sumVectors(aux2,aux4, aux4);    //aux4 = aux4 + aux2
    sumVectors(aux4,aux5, aux5);    //aux5 = aux4 + aux5

    k6 = function(updateCoord(coord, aux5, (double)step, dimension), param);
    std::vector<double> AuxError(dimension);
    for(int i = 0; i < dimension; i++)
    {
        coord[i] = coord[i] + (double)step*(c1*k1[i]+c3*k3[i] + c4*k4[i]+ c5*k5[i]);
        AuxError[i] =  (d1*k1[i]+d3*k3[i] + d4*k4[i]+ d5*k5[i]+d6*k6[i]);
    }
    error = normOf(AuxError);
    if(error>tol)
    {
        coord = coordTemp;
        stepNew = step/2.0;
    }
    else if(error<0.1*tol)
    {
        stepNew = step*2.0;
    }
}
void completeRungeKutta45(std::vector<double> (*function)(std::vector<double>, double), 
                          std::vector<double> InitialCoord,std::vector<std::vector<double>>& rk45Solution, 
                          double param, double step, int dimension,
                          double tol, double maxTime, std::vector<double>& hs)
{
   
    double stepNew, error;
    std::vector<double> tempCoord = InitialCoord;
    rk45Solution.emplace_back(InitialCoord);
  
    double time = 0;
    while(time < maxTime)
    {
        InitialCoord = tempCoord;
        rungeKutta45(function, tempCoord, error ,param, step, dimension, tol, stepNew);
        if(tempCoord==InitialCoord)
        {
            step = stepNew;
            rungeKutta45(function, tempCoord, error ,param, step, dimension, tol, stepNew);
            rk45Solution.emplace_back(tempCoord);
            time += step;
            hs.emplace_back(step);

        }
        else if(stepNew == 2.0*step)
        {
            rk45Solution.emplace_back(tempCoord);
            hs.emplace_back(step);
            time+=step;
            step = stepNew;
            rungeKutta45(function, tempCoord, error ,param, step, dimension, tol, stepNew);
            rk45Solution.emplace_back(tempCoord);
            time+=step;
            hs.emplace_back(step);

        }
        else
        {
            rk45Solution.emplace_back(tempCoord);
            time+=step;
            hs.emplace_back(step);
        }
        // std::cout<<"integrou aqui\n";

    }
}