#ifndef BIFFUCATION_H
#define DIFFURCATION_H
#include<cmath>
#include<iostream>
#include "../include/rungekutta4thSquare.hpp"

std::vector<std::vector<double>> bifurcation(std::vector<double>(*function)(std::vector<double>, double),
                                            double paramRange[2], std::vector<double> initialCond, 
                                            double paramStep, int coordBeingAnalysed  ,int systemDimension);
void printBiffucationToFile(std::vector<std::vector<double>>& matrix, std::string fileName);
void searchMaxMin(std::vector<std::vector<double>>& biffurcation,std::vector<double>& xCoord, std::vector<double>::iterator& paramValue);

#endif