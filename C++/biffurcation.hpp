
#include<cmath>
#include "rungekutta4thSquare.cpp"


std::vector<std::vector<double>> biffurcation(std::vector<double>(*function)(std::vector<double>, double),
                                            double paramRange[2], std::vector<double> initialCond, 
                                            double paramStep, double integrationStep, int dimension);