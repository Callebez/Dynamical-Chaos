#ifndef PRINTING_HPP
#define PRINTING_HPP
#include<vector>
#include<fstream>
#include<iostream>
#include<cmath>

void printLyapunovsToFile(std::vector<std::vector<long double>>& matrix, std::string fileName);
void printMatrix(std::vector<std::vector<double>>& matrix);
void printMatrix(std::vector<std::vector<long double>>& matrix);

void printMatrixToFile(std::vector<std::vector<double>>& matrix,  std::string fileName);
void completeRungeKuttaToFile(std::vector<double> (*function)(std::vector<double>, double), std::vector<double> initialCond,
                                                   double param, double step, int dimension, double time_span[2]);

#endif