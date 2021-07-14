#include"../include/biffurcation.hpp"
#include<vector>

std::vector<double> LyapunovExponents(std::vector<double>(*function)(std::vector<double>,double)
                                      ,std::vector<double>initcoord, double step, int iterations);
//void PrintLastExp(std::vector<double> var1, int times);
void PrintExp(std::vector<double> var1);