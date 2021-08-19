#include "../include/printing.hpp"
#include "../include/rungekutta4thSquare.hpp"

void printLyapunovsToFile(std::vector<std::vector<long double>>& matrix, std::string fileName)
{

    std::vector<std::ofstream> outputs(matrix[0].size());
    std::vector<std::string> files(matrix[0].size());
  
    for(uint i = 0; i < matrix[0].size(); i++)
    {
        files[i]="./outputs/txt/";
        files[i] = files[i] + std::to_string(i) + ".dat";
    }
 
    for(uint i = 0; i < matrix[0].size(); i++)
    {
        outputs[i].open(files[i]);
    }
    for(uint j = 0; j < matrix.size(); j++)
    {
        for(uint i = 0; i < matrix[0].size(); i++)
        {
            outputs[i]<<matrix[j][0]<<" "<< matrix[j][i] <<"\n";
        }
    }
   for(uint i = 0; i < matrix[0].size(); i++)
    {
        outputs[i].close();
    }
 
}
void printMatrix(std::vector<std::vector<double>>& matrix)
{
    for(uint i = 0; i < matrix.size(); i++)
    {
        for(uint j = 0; j < matrix[0].size(); j++)
        {
            std::cout<<matrix[i][j]<<" ";
        }
        std::cout<<"\n";
    }
}
void printMatrix(std::vector<std::vector<long double>>& matrix)
{
    for(uint i = 0; i < matrix.size(); i++)
    {
        for(uint j = 0; j < matrix[0].size(); j++)
        {
            std::cout<<matrix[i][j]<<" ";
        }
        std::cout<<"\n";
    }
}
void printMatrixToFile(std::vector<std::vector<double>>& matrix, std::string fileName)
{
    std::ofstream output;
    std::string file = "./outputs/txt/";
    file = file + fileName;
    output.open(file);
    for(uint i = 0; i < matrix.size(); i++)
    {
        for(uint j = 0; j < matrix[0].size(); j++)
        {
            output<<matrix[i][j]<<" ";
        }
        output<<"\n";
    }
    output.close();
    std::cout<<file;
}
void completeRungeKuttaToFile(std::vector<double> (*function)(std::vector<double>, double), std::vector<double> initialCond,
                                                   double param, double step, int dimension, double time_span[2])
{
    int iterations = (int)(fabs(time_span[1]-time_span[0])/step);
    std::ofstream fileName;
    fileName.open("ouputrk4th.dat");
    //std::vector<double> auxVec = rungeKutta4thSquare(function, initialCond, param, step, dimension);
    for(int j = 0; j < iterations; j++)
    {
        initialCond = rungeKutta4thSquare(function, initialCond, param, step, dimension);
        for(int i = 0; i < dimension; i++)
        {
            fileName << initialCond[i] <<"  ";
        }
        fileName<< std::endl;
    }
    fileName.close();

}