#ifndef LYAPUNOV_H
#define LYAPUNOV_V
#include "../include/LinearAlgebra.hpp"
#include "../include/rungekutta4thSquare.hpp"
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <thread>
std::vector<long double> lyapunovSpectrum(std::vector<double> (*function)(std::vector<double>, double),
                                std::vector<std::vector<double>> (*jacobian)(std::vector<double>&, double, double), 
                                std::vector<double>& initialCond, double tol, double time, 
                                double step, double param);
void laypunovVaringParameter(std::vector<double>(*function)(std::vector<double>, double),
                                            std::vector<std::vector<double>> (*jacobian)(std::vector<double>&, double,double), 
                                            double paramRange[2], std::vector<double> initialCond, double tol,
                                            double paramStep, int coordBeingAnalysed, int systemDimension,
                                            std::vector<std::vector<double>>& lyapunovSpectrumRange);
#endif
