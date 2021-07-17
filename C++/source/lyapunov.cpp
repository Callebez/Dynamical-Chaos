
#include "../include/LyapExp.hpp"
std::vector<double> lorenzJaconian(std::vector<double> coord, std::vector<double> pertubation, double rho)
{
    std::vector<double> jacobian(3,0);

    jacobian[0] = 10.0*pertubation[1]-10.0*pertubation[0];
    jacobian[1] = rho*pertubation[0] - coord[2] * pertubation[0] - pertubation[1] - coord[0]*pertubation[2];
    jacobian[2] = coord[1] * pertubation[0] +coord[0]*pertubation[1] - (8.0/3.0)* pertubation[2];
    return jacobian;
}
std::vector<double> lorenz(std::vector<double>coord, double rho)
{
    std::vector<double> xcoord (3,0);
    xcoord[0] = 10.0*coord[1] - 10.0 * coord[0];
    xcoord[1] = rho*coord[0] - coord[1] - coord[0]*coord[2];
    xcoord[2] = coord[0]*coord[1] -(8.0/3.0)*coord[2];
    return xcoord;
}
double dotProduct(std::vector<double> &v1, std::vector<double> &v2)
{
    double result =  0;
    for(uint i = 0; i < v1.size(); i++)
    {
        result += v1[i]*v2[i];
    }
    return result;
}
double normOf(std::vector<double> &vec)
{
    return sqrt(dotProduct(vec,vec));
}
double normSquare(std::vector<double> &vec)
{
    return dotProduct(vec,vec);
}
void normalize(std::vector<double> &vec)
{
    double norm = normOf(vec);
    for(uint i = 0; i < vec.size(); i ++)
    {
        vec[i]/=norm;
    }
}

double lyapunov(std::vector<double> (*function)(std::vector<double>, double), std::vector<double> initialCond,
                std::vector<double> (*jacobian)(std::vector<double>, std::vector<double>, double), uint dimension ,double step)
{
    uint iterations = (uint)(100.0/step);

    srand(time(NULL));
    std::vector<double> pertubation (dimension,0);
    std::vector<std::vector<double>> coord (iterations+1);
    std::vector<double> xcoord (dimension,0);
    double lyapunov = 0;
    

    for(uint i = 0; i < dimension; i ++)
    {
        pertubation[i] = (double)rand()/RAND_MAX;
    }
    double rho = 28.0;
    coord[0] = rungeKutta4thSquare(function,initialCond, rho,step, dimension);
    pertubation = rungeKutta4thSquarePertubation(jacobian, coord[0], pertubation,rho, step,dimension, xcoord);
    // normalize(xcoord);
    double aux = 0;

    lyapunov += dotProduct(xcoord,pertubation);
    aux += log(normOf(pertubation)); 
    // pertubation = RandVec(pertubation.size());
    // double i = 1.0;
    for(uint i = 0; i < iterations; i++)
    {
        // pertubation = RandVec(pertubation.size());
        normalize(pertubation);

        coord[i+1] = rungeKutta4thSquare(function,coord[i], rho,step, dimension);
        // xcoord = jacobian(initialCond, pertubation, rho);

        pertubation = rungeKutta4thSquarePertubation(jacobian, coord[i], pertubation,rho, step,dimension, xcoord);
        
       
        lyapunov += dotProduct(xcoord,pertubation)/normSquare(pertubation);
        aux += log(normOf(pertubation)); 
        // normalize(xcoord);

        // lyapunov += dotProduct(xcoord,pertubation);
        std::cout<<"novo = "<<lyapunov/(iterations)<<std::endl;
        std::cout<<"velho = "<<aux/(iterations+1)<<"\n"<<std::endl;

        // i++;
    }
    // std::cout<<lyapunov/((double)(iterations));
    return lyapunov/((double)(iterations));

}