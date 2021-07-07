#include "rungekutta4thSquare.hpp"

std::vector<double> rungeKutta4thSquare(std::vector<double> (*function)(std::vector<double>, double), 
                                  std::vector<double> coord, double param, double step, int dimension = 4)
{
    std::vector<double> integratedPoint (dimension,0);
    std::vector<double> k1 (dimension,0);
    std::vector<double> k2 (dimension,0);
    std::vector<double> k3 (dimension,0);
    std::vector<double> k4 (dimension,0);
    
    k1 = function(coord, param);
    k2 = function(updateCoord(coord, k1, step), param);
    k3 = function(updateCoord(coord, k2, step), param);
    k2 = function(updateCoord(coord, k1, step, 1), param);
    for(int i = 0; i <coord.size(); i++)
    {
        integratedPoint[i] = coord[i] + step*(1/6)*(k1[i]+2*k2[i] + 2*k3[i]+ k4[i]);
    }
    return integratedPoint;
      
}
std::vector<double> updateCoord(std::vector<double> coord, std::vector<double> increment, double step, int scale = 2,int dimension = 4)
{
    std::vector<double> updatedCoord (dimension,0);
    for(int i = 0; i < coord.size(); i++)
    {
        updatedCoord[i] = coord[i] + step*increment[i]/scale;
    }
    return updatedCoord;
}