#include"../include/discretelyap.hpp"
#include"../include/LinearAlgebra.hpp"

std::vector<std::vector<double>> discreteSys(std::vector<double> initialcond,double gamma, double tau, int iteration, double k=1)
{
    //Consider (x,y,px,py)
    std::vector < std::vector<double>> auxsys (iteration,std::vector<double>(initialcond.size()));
    auxsys[0] = initialcond;
    for (int i = 1; i < iteration; i++)
    {
        double flu1 = complicate(auxsys[i-1][0],auxsys[i-1][1], gamma, tau);
        double flu2 = complicate(auxsys[i-1][0],auxsys[i-1][1], gamma, tau);
        double forcex = force(auxsys[i - 1][0], auxsys[i - 1][1], k);
        double forcey = force(auxsys[i - 1][1], auxsys[i - 1][0], k);
        //std::cout << flu1 << "   " << flu2 << std::endl;
        auxsys[i][0] = auxsys[i - 1][0] + auxsys[i - 1][2] * tau;
        auxsys[i][1]=auxsys[i-1][1]+auxsys[i-1][3]*tau;
        auxsys[i][2]=(1-gamma*tau)*auxsys[i-1][2]+forcex*tau+flu1;
        auxsys[i][3]=(1-gamma*tau)*auxsys[i-1][3]+forcey*tau+flu2;
    }
    return auxsys;
}
std::vector<long double> discreteLyap(std::vector<double> (*function)(std::vector<double>, double),
                                                   std::vector<std::vector<double>> (*jacobian)(std::vector<double>&, double,double),
                                                   std::vector<std::vector<double>> trajectory, double gamma, double tau)
{
    std::vector<long double> lyapunovExponents (trajectory.size(),0);
    std::vector<double> coord;
    std::vector<std::vector<double>> J = identityMatrix(trajectory[0].size());
    std::vector<std::vector<double>> w = identityMatrix(trajectory[0].size());
    double time = tau * trajectory.size();
    for (uint i = 0; i < trajectory.size(); i++)
    {
        coord = trajectory[i];
        if (i%5000==0)
        {
            std::cout << coord[0] << "   " << coord[1] << "   " << coord[2] << "   " << coord[3] << std::endl;
        }
        J = jacobian(coord, gamma, tau);
        w = matMult(J, w);
        transpostSquare(w);
        gramSchmidt(w);
        for(uint j = 0; j < trajectory[0].size(); j++)
        {
            lyapunovExponents[j] +=  log(normOf(w[j]));
        }
        gramSchmidtNormal(w);
        transpostSquare(w);
    }
    for(uint j = 0; j < trajectory.size(); j++)
    {
        lyapunovExponents[j] = lyapunovExponents[j]/(time);
    }
    return lyapunovExponents;
}
double GaussRand()
{
    unsigned rd = std::chrono::steady_clock::now().time_since_epoch().count();
    //std::random_device rd;
    std::default_random_engine generator (rd);
    std::normal_distribution<double> nd(0, 1);
    double aux = nd(generator);
    return aux;
}

double force(double a, double b, double k)
{
    double result;
    result = -(a + k * a * b * b);
    return result;
}

double complicate(double x, double y, double gamma, double tau)
{
    return GaussRand()*sqrt(2*gamma*tau*laplace(x,y));
}
inline double rho (double x, double y)
{
    double a=0.9122350052;
    double result;
    result = -2 * a + funcrho1(x, y) + funcrho1(y, x) + 2 * (funcrho2(x, y) + funcrho2(y, x));
    return result;
    }
inline double funcrho1(double x, double y)
{
    const double b=0.3345167463;
    const double c=0.0987170706;
    double result;
    result = (b + 6 * c * pow(x, 2)) / (1 + b * (pow(x, 2) + pow(y, 2)) + c * (pow(x, 4) + pow(y, 4)));
    return result;
}
inline double funcrho2(double x, double y)
{
    const double b=0.3345167463;
    const double c=0.0987170706;
    double result;
    result = (b * x + 2 * c * pow(x, 3)) / (1 + b * (pow(x, 2) + pow(y, 2)) + c * (pow(x, 4) + pow(y, 4)));
    result = pow(result, 2);
    return result;
}

double wolfram(double x, double y)
{
    const double a=0.9122350052;
    const double b=0.3345167463;
    const double c=0.0987170706;
    double result;
    result = 8 * (a * pow(1 + b * (pow(x, 2) + pow(y, 2)) + c * (pow(x, 4) + pow(y, 4)), 2) 
           - b * (6 * c * pow(x, 2) * pow(y, 2) + 1)
           + c * (pow(x, 2) + pow(y, 2))
          * (c * (pow(x, 4) - 4 * pow(x, 2) * pow(y, 2) + pow(y, 4)) - 3));
    result /= pow(b * (pow(x, 2) + pow(y, 2)) + c * (pow(x, 4) + pow(y, 4)) + 1, 2);
    result *= -1;
    //std::cout << result << std::endl;
    return result;
}


double laplace(double x, double y)
{
    const double a=0.9122350052;
    const double b=0.3345167463;
    const double c=0.0987170706;
    double result = 0;
    result += -4 * a + 2 * ((2 * b + 12 * c * x) / u(x, y) + pow((2 * b * x + 4 * c * pow(x, 3)) / u(x, y), 2));
    result += -4 * a + 2 * ((2 * b + 12 * c * y) / u(x, y) + pow((2 * b * y + 4 * c * pow(y, 3)) / u(x, y), 2));
    return -result/4;
}
double u(double x, double y)
{
    const double b=0.3345167463;
    const double c=0.0987170706;
    double result;
    result = 1 + b * (pow(x, 2) + pow(y, 2)) + c * (pow(x, 4) + pow(y, 4));
    return result;
}