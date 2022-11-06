#include "numericc.h"

double f1(double x)
{
    return exp(x) - 2 * x * cos(x) - 3;
}
double df1(double x)
{
    return exp(x) - 2 * cos(x) + 2 * x * sin(x);
}

double f2(double x)
{
    return pow(x, 2) + sin(x) + exp(x) - 2;
}

double df2(double x)
{
    return 2 * x + cos(x) + exp(x);
}

double g(double x)
{
    return -sqrt((2.0 / 3.0) * exp(x));

}


double Lag1(double x)
{
    return 1 + sin(x*3.14);
}
double Lag2(double x)
{
    return 1 / (1 + pow(x, 2));
}

using RFS = Numericc::Solutions::RootFindingSolution;
using LSS = Numericc::Solutions::LinearSystemSolution<double>;
using PIS = Numericc::Solutions::PowerIterationSolution<double>;

using integerMat = Numericc::Super::Matrix<int>;

using namespace Numericc::RootFinding;
using namespace Numericc::Interpolators::Newton;
using namespace Numericc::Interpolators;
using namespace Numericc::LinearAlgebra;

int main()
{
    //RFS fixed = FixedPoint(g, -0.5, 100, 5);
    //RFS bisection = Bisection(f1, 0, 2, 5, 5);
    //RFS newton = NewtonRaphson(f1, df1, 1, 10, 5);


    //std::vector<double> x3 = { 45, 50, 55, 60 };
    //std::vector<double> y3 = { 0.7071 , 0.7660,0.8192,0.8660 };

    //std::cout << NewtonGregoryBD(52, x3, y3) << std::endl;
    //std::cout << NewtonGregoryFD(52, x3, y3) << std::endl;
    //std::cout << Lagrange(52, x3, y3) << std::endl;


    integerMat A = { {4, 1, 2, 1},{1, 7, 1, 0},{2, 1, 4, -1},{1, 0, -1, 3} };
    integerMat B = { -8,-20,-2,4 };

    LSS linear = LinearSystemSolver(A, B);
    //PIS pis = PowerIteration(A, 100);
    std::cout << linear.A << " " << linear.B << " " << linear.X;
    std::cout << A << B;

    std::cout << A.Transpose();

    //std::cout << pis.eigenvalue<< " " << pis.eigenvector;

    return 0;
   
}