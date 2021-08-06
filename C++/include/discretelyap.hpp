#include<vector>
#include<cmath>
#include<ctime>
#include<iostream>
#include<random>
#include<chrono>
std::vector<std::vector<double>> discreteLyap(std::vector<double> initialcond, double gamma, double tau, int iteration, double k);
double GaussRand();
double force(double a, double b, double k);
double complicate(double x, double y, double gamma,double tau  );
inline double rho (double x, double y);
inline double funcrho1 (double x, double y);
inline double funcrho2 (double x, double y);
