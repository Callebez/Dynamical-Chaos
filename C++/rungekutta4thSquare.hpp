#include<vector>

std::vector<double> updateCoord(std::vector<double> coord, std::vector<double> increment, double step, int scale = 2,int dimension = 4);
std::vector<double> rungeKutta4thSquare(std::vector<double> (*function)(std::vector<double>, double), 
                                  std::vector<double> coord, double param, double step, int dimension = 4);
