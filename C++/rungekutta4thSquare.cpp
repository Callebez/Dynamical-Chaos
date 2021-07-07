#include "rungekutta4thSquare.hpp"
#include <fstream>

// void printVector(std::vector<double> V, std::ofstream &ouput)
// {
//     for(int i = 0; i < V.size(); i++)
//     {
//         ouput << V[i] << std::endl;
//     }
// }
std::vector<double> rungeKutta4thSquare(std::vector<double> (*function)(std::vector<double>, double), 
                                  std::vector<double> coord, double param, double step, 
                                  int dimension )
{
    std::vector<double> k1 (dimension);
    std::vector<double> k2 (dimension);
    std::vector<double> k3 (dimension);
    std::vector<double> k4 (dimension);
    
    k1 = function(coord, param);
    k2 = function(updateCoord(coord, k1, step, 2,dimension), param);
    k3 = function(updateCoord(coord, k2, step, 2,dimension), param);
    k4 = function(updateCoord(coord, k1, step, 1,dimension), param);
    for(int i = 0; i < dimension; i++)
    {
        coord[i] = coord[i] + step/6.0*(k1[i]+2*k2[i] + 2*k3[i]+ k4[i]);
    }
    
    // for(int i = 0; i < dimension; i++)
    // {
    //     std::cout << "Coord: " << coord[i]<<". k1 "<< k1[i]<< ", k2: "<<k2[i]<< ", k3: "<< k3[i]<< ", k4: "<< k4[i]<< std::endl;
    // }
    return coord;
      
}
std::vector<double> updateCoord(std::vector<double> coord, std::vector<double> increment, double step, int scale,int dimension)
{
    std::vector<double> updatedCoord (dimension);
    for(int i = 0; i < dimension; i++)
    {
        updatedCoord[i] = coord[i] + step*increment[i]/scale;
    }
    return updatedCoord;
}