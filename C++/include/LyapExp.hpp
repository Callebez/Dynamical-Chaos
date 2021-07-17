#ifndef LYAPUNOV_H
#define LYAPUNOV_V
#include "../include/rungekutta4thSquare.hpp"
#include <ctime>
#include <cstdlib>
#include<iostream>
std::vector<double> lorenzJaconian(std::vector<double> coord, std::vector<double> pertubation, double rho);
std::vector<double> lorenz(std::vector<double>coord, double rho);
double dotProduct(std::vector<double> &v1, std::vector<double> &v2);
double normOf(std::vector<double> &vec);
void normalize(std::vector<double> &vec);
double lyapunov(std::vector<double> (*function)(std::vector<double>, double), std::vector<double> initialCond,
                std::vector<double> (*jacobian)(std::vector<double>, std::vector<double>, double), uint dimension ,double step);
std::vector<long double> RandVec(uint dimention);

#endif