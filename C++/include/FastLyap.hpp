#include<vector>

std::vector<double> RandVec(int dimention);
std::vector<double> normalize (std::vector<double> vec);
long double FastLyapExp (std::vector<double> (*system)(std::vector<double>, double),
                     std::vector<std::vector<double>> (*jacobian)(std::vector<double>, double),
                     double parameter,std::vector<double> initialcondition, double step, int iterations);
long double dotproduct (std::vector<double> a , std::vector<double> b);
std::vector<double> MatrixVector (std::vector<std::vector<double>> mat, std::vector<double> vec);
long double VecLenght (std::vector<double> vec);