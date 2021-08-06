#include"../include/discretelyap.hpp"

std::vector<std::vector<double>> discreteLyap(std::vector<double> initialcond,double gamma, double tau, int iteration, double k=1)
{
    //Consider (x,y,px,py)
    std::vector < std::vector<double>> auxsys (iteration,std::vector<double>(initialcond.size()));
    auxsys[0] = initialcond;
    for (int i = 1; i < iteration; i++)
    {
        double flu = complicate(auxsys[i-1][0],auxsys[i-1][1], gamma, tau);
        double forcex = force(auxsys[i - 1][0], auxsys[i - 1][1], k);
        double forcey = force(auxsys[i - 1][1], auxsys[i - 1][0], k);

        auxsys[i][0]=auxsys[i-1][0]+auxsys[i-1][2]*tau;
        auxsys[i][1]=auxsys[i-1][1]+auxsys[i-1][3]*tau;
        auxsys[i][2]=(1-gamma*tau)*auxsys[i-1][2]-forcex*tau+flu*GaussRand();
        auxsys[i][3]=(1-gamma*tau)*auxsys[i-1][3]-forcey*tau+flu*GaussRand();
    }
    return auxsys;
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

    return sqrt(2*gamma*tau*rho(x,y));
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