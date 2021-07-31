#ifndef LORENZEQUATIONS_HPP
#define LORENZEQUATIONS_HPP
#include <vector>
#include "../include/LinearAlgebra.hpp"
std::vector<std::vector<double>> lorenzJacobian(std::vector<double>& coord, double rho, double step);
std::vector<double> lorenz(std::vector<double> coord, double rho);

#endif
