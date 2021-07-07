#include <vector>
#include <stdlib.h>
#include <time.h>
#include <cmath>

std::vector<double> classicalPendulum(std::vector<double> coord, double k);
std::vector<double> quantumPenduli(std::vector<double> coord, double gamma);
double flutuation(double x, double y, double gamma);