#include <iostream> //Permite a entrada e saída de valores armazenados, ou a ser armazenados no próprio programa
#include <fstream> //Permite a entrada e saída de valores em documentos externos ao programa
#include <string> //Permite ações com valores de caractere
#include <vector> //Permite a implementação de vetores
#include <cmath> //Libraría numérica para operações matemáticas
#include <stdio.h>
#include <iomanip>
#include <conio.h> //Console Input-Output library
#include <ctype.h> //Libraría para tipos de funções
using namespace std; //Determina qual libraria de fontes utilizar

//Arquivo de funções para programa 5.1
void euler_exp(double (*der)(double),double xini,double xfin,double y,double inc,vector<double> &xvec,vector<double> &yvec) //Solução numérica por método de Euler
{
    double x=xini;
    xvec.push_back(x);
    yvec.push_back(y);
    while(x<=xfin)
    {
        y+=der(y)*inc;
        x+=inc;
        xvec.push_back(x);
        yvec.push_back(y);
    }
}
void euler_imp(double (*der)(double),double xini,double xfin,double y,double inc,vector<double> &xvec,vector<double> &yvec) //Solução numérica por método de Euler
{
    double x=xini;
    double sol[2]={y,0};
    while(x<=xfin)
    {
        sol[1]=sol[0]+der(sol[1])*inc;
        sol[0]=sol[1];
        x+=inc;
        xvec.push_back(x);
        yvec.push_back(sol[1]);
    }
}
void data1(int i,vector<double> &x, vector<double> &y)
{
    string val=to_string(i);
    ofstream dat1;
    dat1.open("Data/exp/valores"+val+".dat");
    for(int i=0;i<x.size()-1;i++)
    {
        dat1<<x[i]<<"    "<<y[i]<<endl;
    }
    dat1<<endl;
    dat1.close();  
}
void data2(int i,vector<double> &x, vector<double> &y)
{
    string val=to_string(i);
    ofstream dat2;
    dat2.open("Data/imp/valores"+val+".dat");
    for(int i=0;i<x.size()-1;i++)
    {
        dat2<<x[i]<<"    "<<y[i]<<endl;
    }
    dat2<<endl;
    dat2.close();
}
void plot(int i)
{
    FILE *gnupipe = _popen("gnuplot -persist", "w");
    if(gnupipe)
    {
        fprintf(gnupipe,"set terminal pdf size 10,7 \n");
        fprintf(gnupipe,"set lmargin 12  \n");
        if(i==0)
        {
            fprintf(gnupipe, "cd \"Data/exp\"\n");
            fprintf(gnupipe,"set output \'../../Graph/graph1.pdf\'\n");
            fprintf(gnupipe,"set title \"Practice 5.1.1\" font \"Times New Roman,40\" \n");
            fprintf(gnupipe,"set label \'Método explícito\' at screen 0.525,0.9 center font \"Times New Roman,30\" \n");
        }
        else
        {
            fprintf(gnupipe, "cd \"Data/imp\"\n");
            fprintf(gnupipe,"set output \'../../Graph/graph2.pdf\'\n");
            fprintf(gnupipe,"set title \"Practice 5.1.2\"font \"Times New Roman,40\" \n");
            fprintf(gnupipe,"set label \'Método implícito\' at screen 0.525,0.9 center font \"Times New Roman,30\" \n");
        }
        fprintf(gnupipe,"set tics font \",20\"\n");
        fprintf(gnupipe,"set xtics 0,5\n");
        fprintf(gnupipe,"set key font \',20\'\n");
        fprintf(gnupipe,"set xlabel \"X\" font\"Times New Roman,30\" offset 0,-1\n");
        fprintf(gnupipe,"set ylabel \"f(x)\" font\"Times New Roman,30\" offset -1,0\n");
        fprintf(gnupipe,"set grid lw 2 \n");
        fprintf(gnupipe,"plot \"valores1.dat\" with line title \"Incremento de 5\",\"valores2.dat\" with line  title \"Incremento de 0.5\",\"valores3.dat\" with line  title \"Incremento de 0.05\",\"valores4.dat\" with line  title \"Incremento de 0.005\",1000*exp(-0.1*x) title \"Solução analítica\" \n      ");        

    }
}
