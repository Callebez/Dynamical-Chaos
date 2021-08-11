#include<vector>
#include<cmath>
#include<ctime>
#include<iostream>
#include<random>
#include<chrono>
std::vector<std::vector<double>> discreteSys(std::vector<double> initialcond, double gamma, double tau, int iteration, double k);
std::vector<long double> discreteLyap(std::vector<double> (*function)(std::vector<double>, double),
                                      std::vector<std::vector<double>> (*jacobian)(std::vector<double> &, double, double),
                                      std::vector<std::vector<double>> trajectory, double gamma, double tau);
double GaussRand();
double force(double a, double b, double k);
double complicate(double x, double y, double gamma,double tau  );
inline double rho (double x, double y);
inline double funcrho1 (double x, double y);
inline double funcrho2 (double x, double y);
double wolfram(double x, double y);
double laplace(double x, double y);
double u(double x, double y);