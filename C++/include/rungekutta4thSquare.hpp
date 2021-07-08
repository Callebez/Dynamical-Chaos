#include <fstream>
#include <string>
# include <cstdlib>

#ifndef RUNGEKUTTA4TH_H
#define RUNGEKUTTA4TH_H
#include<vector>

std::vector<double> updateCoord(std::vector<double> coord, std::vector<double> increment, double step, int scale = 2,int dimension = 4);
std::vector<double> rungeKutta4thSquare(std::vector<double> (*function)(std::vector<double>, double), 
                                  std::vector<double> coord, double param, double step, int dimension);
void completeRungeKuttaToFile(std::vector<double> (*function)(std::vector<double>, double), std::vector<double> initialCond,
                                                   double param, double step, int dimension, double time_span[2]);

#endif