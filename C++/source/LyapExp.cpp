#include"../include/biffurcation.hpp"
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
        norm[0]+=pow(coord[l][0],2);
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
            norm[i]+=pow(coord[l][i],2);
            norm[i]=sqrt(norm[i]);
            coord[l][i]/=norm[i];
        }
            //Cálculo de expoente de Lypunov
        for(int j=0;j<NumEq;j++)
        {
            exponents[j]=+log(norm[j])/log(2.0);
        }
        std::cout<<coord[l][0]<<"    "<<exponents[0]/coord[l][0]
                              <<"    "<<exponents[1]/coord[l][0]
                              <<"    "<<exponents[2]/coord[l][0]
                              <<"    "<<exponents[3]/coord[l][0]<<std::endl;
    }
    /*for(int l=0;l<10;l++)
    {
        for(int i=0;i<NumEq;i++)
        {
            std::cout<<coord[l][i]<<"   ";
        }   
        std::cout<<std::endl;
    }*/
    return exponents;
}