#include"../include/biffurcation.hpp"
#include<vector>

std::vector<double> LyapunovExponents(std::vector<double>(*function)(std::vector<double>,double)
                                      ,std::vector<double>initcoord, double step, int iterations);