#include <vector>
#include <cstdlib>
#include<ctime>
#include <cmath>

#ifndef PENDULI_H
#define PENDULI_H

std::vector<double> classicalPendulum(std::vector<double> coord, double k);
std::vector<double> quantumPendulum(std::vector<double> coord, double gamma);
double flutuation(double x, double y, double gamma);
void waveFunction(double x, double y, double& rho);

std::vector<std::vector<double>> classicalPendulumJacobian(std::vector<double>& coord, double k, double step);
#endif