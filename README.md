# numericc

A small C++ header-only library of numerical methods for linear algebra, root-finding, interpolators and more.

## Examples

### root-finding

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
