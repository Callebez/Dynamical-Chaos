#ifndef LYAPUNOV_H
#define LYAPUNOV_V
#include "../include/LinearAlgebra.hpp"
#include "../include/rungekutta4thSquare.hpp"
#include <ctime>
#include <cstdlib>
#include<iostream>

std::vector<long double> lyapunovSpectrum(std::vector<double> (*function)(std::vector<double>, double),
                                std::vector<std::vector<double>> (*jacobian)(std::vector<double>&,double), 
                                std::vector<double>& initalCond, double step, double param);
#endif