#include <fstream>



#ifndef RUNGEKUTTA4TH_H
#define RUNGEKUTTA4TH_H
#include<vector>
#include<cmath>

std::vector<double> updateCoord(std::vector<double>coord, std::vector<double> increment, double step, double scale,int dimension);
std::vector<double> rungeKutta4thSquare(std::vector<double> (*function)(std::vector<double>, double), 
                                  std::vector<double>& coord, double param, double step, 
                                  int dimension );
std::vector<double> rungeKutta4thSquarePertubation(std::vector<double> (*function)(std::vector<double>,std::vector<double>, double), 
                                  std::vector<double> coord,std::vector<double> pertubation, double param, double step, 
                                  int dimension, std::vector<double>& functionAval );
#endif