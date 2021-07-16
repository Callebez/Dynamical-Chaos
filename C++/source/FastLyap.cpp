#include<cmath>
#include<vector>
#include <time.h>
#include<iostream>
#include"../include/rungekutta4thSquare.hpp"

std::vector<double> RandVec(int dimention)
{
    std::vector<double> z;
    srand(time(NULL));
    for(int i=0;i<dimention;i++)
    {
        z.push_back(((long double)rand()/RAND_MAX)*pow(-1,rand()%2));
    }
    std::cout<<rand()<<std::endl<<rand()<<std::endl<<rand()<<std::endl;
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
long double LyapExp (std::vector<double> (*function)(std::vector<double>, double),std::vector<double> initialcondition, double step)
{
    long double exponent;
    

    return exponent;
}
long double dotproduct (std::vector<double> a , std::vector<double> b)
{
    long double result = 0.0;
    if(a.size()==b.size())
    {
        for(int i=0;i<a.size();i++)
        {
            result+=a[i]*b[i];
        }
    }
    return result;
}

