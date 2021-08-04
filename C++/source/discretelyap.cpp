#include"../include/discretelyap.hpp"

std::vector<std::vector<double>> discreteLyap(std::vector<double> initialcond,double gamma, double tau, int iteration, double k=1)
{
    //Consider (x,y,px,py)
    std::vector < std::vector<double>> auxsys (iteration,std::vector<double>(initialcond.size()));
    auxsys[0] = initialcond;
    for (int i = 1; i < iteration; i++)
    {
        double f = flutuation(auxsys[i-1][0],auxsys[i-1][1], gamma, tau);
/*
Colocar na função flutuation as características necessárias 
*/
        auxsys[i][0]=auxsys[i-1][0]+auxsys[i-1][2]*tau;
        auxsys[i][1]=auxsys[i-1][1]+auxsys[i-1][3]*tau;
        auxsys[i][2]=auxsys[i-1][2]+(1-gamma*tau)*auxsys[i-1][2]+force(auxsys[i-1][0],auxsys[i-1][1],k)*tau+f;
        auxsys[i][3]=auxsys[i-1][3]+(1-gamma*tau)*auxsys[i-1][3]+force(auxsys[i-1][1],auxsys[i-1][0],k)*tau+f;
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

static double force(double a, double b, double k)
{
    double result;
    result = -(a + k * a * b * b);
    return result;
}

static double flutuation(double x, double y, double gamma, double tau)
{
    
    return  GaussRand()*sqrt(2*(gamma)*tau*((0.62832912000*(pow(x,2.0) + pow(y,2.0)) +
                        0.01205153100* pow(x,2.0) * pow(y,2.0)  +
                        0.56437351570*(pow(x,4.0) + pow(y,4.0)) +
                        0.13998728990*(pow(x,6.0) + pow(y,6.0)) +
                        0.06202680930*(pow(x,2.0)* pow(y,4.0) +  pow(y,2.0)*pow(x,4.0)) +
                        0.03555913955* pow(x,4.0) *  pow(y,4.0)  +
                        0.01777956970*( pow(x,8.0) +  pow(y,8.0)) + 1.155436517)/
                        pow((1.0 + 0.3345167463*(pow(x,2.0) +  pow(y,2.0)) +
                        0.09871707060*(pow(x,4.0) +  pow(y,4.0))),2.0)));
}


