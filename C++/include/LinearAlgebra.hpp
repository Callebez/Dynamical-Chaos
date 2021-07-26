#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP
#include <vector>
#include <fstream>
#include <cmath>
#include <string>

//Very basic operations
inline long double dotProduct(std::vector<double> &v1, std::vector<double> &v2);
long double normOf(std::vector<double> &vec);
inline double normSquare(std::vector<double> &vec);
inline void normalize(std::vector<double> &vec);

//Operations on matrices
void gramSchmidt(std::vector<std::vector<double>>& matrix);
void gramSchmidtNormal(std::vector<std::vector<double>>& matrix); 
void transpostSquare(std::vector<std::vector<double>>& matrix);

//auxiliary operations on matrices
std::vector<double> projectionIntoU(std::vector<double> vectorV,std::vector<double> vectorU);
std::vector<double>  sumOfProjections(uint position, std::vector<std::vector<double>>& matrix);

std::vector<std::vector<double>> matMult(std::vector<std::vector<double>>& leftMatrix, std::vector<std::vector<double>>& rigthMatrix);
std::vector<std::vector<double>> identityMatrix(uint order);


#endif