# numericc

A small C++ header-only library of numerical methods for linear algebra, root-finding, interpolators and ODE solving.

### Contents:
* [Examples](#examples)
    * [Linear Algebra](#linear-algebra)
    * [Root finding](#root-finding)
    * [ODE solving](#ode-solving)
    * [Interpolation](#interpolation)
* [Installation and Dependencies](#installation-and-dependencies)
* [Author and License](#author)

## Examples

### Linear Algebra

```cpp
using integerMat = Numericc::Super::Matrix<int>;

using namespace Numericc::LinearAlgebra;

integerMat A = { {4, 1, 2, 1}, {1, 7, 1, 0}, {2, 1, 4, -1}, {1, 0, -1, 3} };
integerMat C = { {1, 8, -2, 8}, {10, 1, -1, 0}, {-4, -6, 16, 0}, {10, -9, -7, 1} };
integerMat B = { -8, -20, -2, 4}; // 1D matrix representation

/* Matrix manipulations */

auto A_trans = A.Transpose();
auto A_det = A.Determinant();
auto A_inverse = A.Inverse();
auto A_adjoint = A.Adjoint();
auto A_minmax = A.MinMax(); // struct {T min,max};

auto D = A * C // matrix multiplication
auto E = D / 5;

/* Linear Algebra */

/* struct { Matrix<T> A, Matrix<T> B // components of the resulting A|B triangular augmented matrix.                                                                             , Matrix<T> X              // solution to Ax = B } */
using LSS = Numericc::Solutions::LinearSystemSolution<double>;

/* struct  { T eigenvalue                   // largest eigen value
           , std::vector<T> eigenvector     // corresponding eigenvector } */                                                                      
using PIS = Numericc::Solutions::PowerIterationSolution<double>; // struct

LSS linear = LinearSystemSolver(A, B); // solves linear system Ax = B .

auto correct = A*linear.X == B; // true
std::cout << linear.A << linear.B << linear.A << std::endl; // overloaded << for matrices .

size_t max_iters = 100;
PIS pis = PowerIteration(A, max_iters); // computes largest eigenvalue of matrix and corresponding eigenvector.

std::cout << pis.eigenvalue << " " << pis.eigenvector << std::endl;
}
```
### Root finding

```cpp
    double f(double x)
    {
        return exp(x) - 2 * x * cos(x) - 3;
    }
    double df(double x)
    {
        return exp(x) - 2 * cos(x) + 2 * x * sin(x);
    }

    using RFS = Numericc::Solutions::RootFindingSolution; // struct{double root, size_t iterations}

    using namespace Numericc::RootFinding;

    std::vector<double> bracket = {0,2};
    auto a = bracket[0];
    auto b = bracket[1];
    size_t max_iterations = 50;
    size_t accuracy = 5; // Scarborough criterion
    RFS bisection = Bisection(f, a, b, max_iterations, accuracy);

    double guess = 1.0;
    RFS newton = NewtonRaphson(f, df, guess, max_iterations, accuracy);
```

### ODE solving

```cpp
double RHS(double x, double y)
{
    return x + y + x * y;
}

/* struct { double X, Y    // y (X) at X
            double Xs, Ys  // solution vectors up to X} */
using ODES = Numericc::Solutions::OdeSolution;

using namespace Numericc::ODE;
using namespace Numericc::ODE::RK;

double x0 = 0.0;
double y0 = 1.0;
double step = 0.025;
double x = 0.54;

ODES euler = Euler(x0, y0, step, x, RHS);
ODES rk4 = RK4(x0, y0, step, x, RHS);
```

### Interpolation

```cpp
double f(double x)
{
    return 1 / (1 + pow(x, 2));
}

using namespace Numericc::Interpolators::Newton;
using namespace Numericc::Interpolators;

std::vector<double> xv = { 45, 50, 55, 60 };
std::vector<double> yv = { 0.7071, 0.7660, 0.8192, 0.8660 };

double x = 52;
double lagrange =  Lagrange(x, xv, yv);
double newtonBD =  NewtonGregoryBD(x, xv, yv);
double newtonFD =  NewtonGregoryFD(x, xv, yv)
```

## Installation and Dependencies
numericc is a header-only library. Simply add the header files to your project using
```cpp
#include "include/numericc.h"
```
The only dependency is a C++11 compatible compiler.

## Author

Evangelos Smyrniotis

## License

MIT License
