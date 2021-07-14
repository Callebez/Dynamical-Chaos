#include"../include/biffurcation.hpp"
#include"../include/LyapExp.hpp"
#include<cmath>
#include<iostream>

std::vector<double> LyapunovExponents (std::vector<double>(*function)(std::vector<double>,double)
                                      ,std::vector<double>initcoord, double step, int iterations)
{
    int NumEq=initcoord.size(); 
    std::vector<double> exponents(NumEq,0);
    std::vector<double> norm(NumEq,0);
    std::vector<double> GramCoef(NumEq,0);
    std::vector<std::vector<double>> coord(iterations);
    coord[0]=initcoord;
    for(int l=1;l<iterations;l++)
    {
        coord[l]=rungeKutta4thSquare(function, coord[l-1], 1e-3, step, NumEq);
        //Normalização do primeiro vetor
       /* for(int i=0;i<NumEq;i++)
        {
            norm[0]+=pow(coord[l][NumEq*i+1],2);
        }
        norm[0]=sqrt(norm[0]);
        for(int i=0;i<NumEq;i++)
        {
            coord[l][NumEq*i+1]/=norm[0];        
        }*/
        for(int i=0;i<NumEq;i++)
        {
            norm[0]+=pow(coord[l][i],2);
        }
        norm[0]=sqrt(norm[0]);
        coord[l][0]/=norm[0];        
        //Normalização dos demais vetores
        for(int i=1;i<NumEq;i++)
        {
            for(int j=0;j<i-1;j++)
            {
                    GramCoef[j]+=coord[l][i]*coord[l][j];
            }
            for(int k=0;k<i-1;k++)
            {
                coord[l][i]-=GramCoef[k]*coord[l][i];
            }
            for(int k=0;k<NumEq;k++)
            {
                norm[i]+=pow(coord[l][i],2);
            }
            norm[i]=sqrt(norm[i]);
            coord[l][i]/=norm[i];
        }
            //Cálculo de expoente de Lypunov
        for(int j=0;j<NumEq;j++)
        {
            exponents[j]=+log(norm[j])/log(2.0);
            exponents[j]/=coord[l][0];
        }
        PrintExp(exponents);
    }

    std::cout<<std::endl;
    //PrintLastExp(exponents,20);
    return exponents;
}
/*
void PrintLastExp(std::vector<double> var1, int times)
{
    for(int j=0;j>=times;j--)
    {
        for(int i =0;i<(int)var1.size();i++)
        {
            std::cout<<var1[i]<<"    ";
        } 
        std::cout<<std::endl;
    }
}
*/
void PrintExp(std::vector<double> var1)
{
    for(int i=0;i<(int)var1.size();i++)
        {
            std::cout<<var1[i]<<"    ";
        } 
        std::cout<<std::endl;
}