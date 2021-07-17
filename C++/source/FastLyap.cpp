#include<cmath>
#include<vector>
#include <time.h>
#include<iostream>
#include"../include/rungekutta4thSquare.hpp"
#include"../include/systems.hpp"
#include"../include/FastLyap.hpp"

std::vector<double> RandVec(int dimention)
{
    std::vector<double> z;
    srand(time(NULL));
    for(int i=0;i<dimention;i++)
    {
        z.push_back(((long double)rand()/RAND_MAX)*pow(-1,rand()%2));
    }
    /*for(int i=0;i<dimention;i++)
    {
        std::cout<<z[i]<<"   ";
    }*/
    return z;
}
std::vector<double> normalize (std::vector<double> vec)
{
    double mod=0;
    for(uint i=0;i<(uint)vec.size();i++)
    {
        mod+=pow(vec[i],2);
    }
    mod=sqrt(mod);
    for(uint i=0;i<(uint)vec.size();i++)
    {
        vec[i]/=mod;
    }
    return vec;
}
long double FastLyapExp (std::vector<double> (*system)(std::vector<double>, double),
                     std::vector<std::vector<double>> (*jacobian)(std::vector<double>, double),
                     double parameter,std::vector<double> initialcondition, double step, int iterations)
{
    long double exponent;
    std::ofstream dados("./DAT/Fast/LyapunovExponents.dat");
    std::vector<std::vector<double>> coord (iterations);
    coord[0]=initialcondition;
    std::vector<double>Z ((int)initialcondition.size());
    std::vector<double>UZ ((int)initialcondition.size());
    for(int i=0;i<iterations-1;i++)
    {
        Z=normalize(RandVec((int)Z.size()));
        UZ=MatrixVector(LorenzJacobian(coord[i],parameter),Z);
        exponent=dotproduct(UZ,Z)/pow(VecLenght(Z),2);
        dados<<i<<"   "<<exponent<<std::endl;
        coord[i+1]=rungeKutta4thSquare(system,coord[i],parameter,step,(int)coord[0].size());
    }
    return exponent;
}
long double dotproduct (std::vector<double> a , std::vector<double> b)
{
    long double result = 0.0;
    if(a.size()==b.size())
    {
        for(int i=0;i<(int)a.size();i++)
        {
            result+=a[i]*b[i];
        }
    }
    return result;
}
std::vector<double> MatrixVector (std::vector<std::vector<double>> mat, std::vector<double> vec)
{
    std::vector<double> result;
    double aux;
    if(mat[0].size()==vec.size())
    {
        for(int i=0;i<(int)mat.size();i++)
        {
            aux=0;
            for(int j=0;j<(int)vec.size();j++)
            {
                aux+=mat[i][j]*vec[j];
            }
            result.push_back(aux);
        }
    }
    return result;
}
long double VecLenght (std::vector<double> vec)
{
    long double lenght=0;
    for (int i=0;i<(int)vec.size();i++)
    {
        lenght+=pow(vec[i],2);
    }
    lenght=sqrt(lenght);
    return lenght;
}




