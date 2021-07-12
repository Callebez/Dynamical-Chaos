#include"../include/LyapExp.hpp"
#include<cmath>
#include<iostream>

std::vector<double> LyapunovExponents (std::vector<double>coord)
{
    int NumEq=4; 
    std::vector<double> exponents(NumEq,0);
    std::vector<double> norm(NumEq,0);
    std::vector<double> GramCoef(NumEq,0);
    //Normalização do primeiro vetor
    for(int i=0;i<NumEq;i++)
    {
        norm[0]+=pow(coord[NumEq*i+1],2);
    }
    norm[0]=sqrt(norm[0]);
    for(int i=0;i<NumEq;i++)
    {
        coord[NumEq*i+1]/=norm[0];        
    }
    //Normalização dos demais vetores
    for(int i=1;i<NumEq;i++)
    {
        for(int j=0;j<i-1;j++)
        {
            for(int k=0;k<NumEq;k++)
            {
                GramCoef[j]+=coord[NumEq*k+i]*coord[NumEq*k+j];
            }
        }
        for(int j=0;j<NumEq;j++)
        {
            for(int k=0;k<i-1;k++)
            {
                coord[NumEq*j+i]-=GramCoef[k]*coord[NumEq*j+i];
            }
        }
        for(int j=0;j<NumEq;j++)
        {
            norm[i]+=pow(coord[NumEq*j+i],2);
        }
        norm[i]=sqrt(norm[i]);
        for(int j=0;j<NumEq;j++)
        {
            coord[NumEq*j+i]*=norm[i];
        }
        //Cálculo de expoente de Lypunov
        for(int j=0;j<NumEq;j++)
        {
            exponents[j]=+log(norm[j])/log(2.0);
        }
        
        
        /*Mensagem para Callebe:
        No programa do cálculo dos expoentes de Lyapunov em Fortran
        Eles dividem pelo logaritmo neperiano de 2.
        Mas, ao que eu entendo esse valor deveria ser calculado
        a cada vez que este programa for rodado.
        Não sei se minha dúvida é boa, ou é somente
        meu mal entendimento.*/

        //Falta realizar o Loop com label 100 no programa original de Fortran
        //Junto com a invocação da função de integração  

        //std::cout<<exponents[];


    }   








    return exponents;
}