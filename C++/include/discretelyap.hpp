#include<vector>
#include<cmath>
#include<ctime>
#include<iostream>
#include<random>
#include<chrono>
std::vector<std::vector<double>> discreteLyap(std::vector<double> initialcond, double gamma, double tau, int iteration, double k);
double GaussRand();
static double force(double a, double b, double k);
static double flutuation(double x, double y, double gamma,double tau  );
