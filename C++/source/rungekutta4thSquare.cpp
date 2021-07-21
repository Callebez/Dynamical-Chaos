#include "../include/rungekutta4thSquare.hpp"


std::vector<double> rungeKutta4thSquare(std::vector<double> (*function)(std::vector<double>, double), 
                                  std::vector<double> &coord, double param, double step, 
                                  int dimension )
{
    std::vector<double> k1 (dimension);
    std::vector<double> k2 (dimension);
    std::vector<double> k3 (dimension);
    std::vector<double> k4 (dimension);
    
    k1 = function(coord, param);
    k2 = function(updateCoord(coord, k1, (double)step, 2.0,dimension), param);
    k3 = function(updateCoord(coord, k2, (double)step, 2.0,dimension), param);
    k4 = function(updateCoord(coord, k3, (double)step, 1.0,dimension), param);
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

void completeRungeKuttaToFile(std::vector<double> (*function)(std::vector<double>, double), std::vector<double> initialCond,
                                                   double param, double step, int dimension, double time_span[2])
{
    int iterations = (int)(fabs(time_span[1]-time_span[0])/step);
    std::ofstream fileName;
    fileName.open("ouputrk4th.dat");
    //std::vector<double> auxVec = rungeKutta4thSquare(function, initialCond, param, step, dimension);
    for(int j = 0; j < dimension; j++)
    {
        initialCond = rungeKutta4thSquare(function, initialCond, param, step, dimension);
        for(int i = 0; i < iterations; i++)
        {
            fileName << initialCond[j] <<"  ";
        }
        fileName<< std::endl;
    }
    fileName.close();

}