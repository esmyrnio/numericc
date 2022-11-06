# numericc

A small C++ header-only library of numerical methods for linear algebra, root-finding, interpolators and more.

# Examples

```cpp
    using RFS = Numericc::Solutions::RootFindingSolution;
    using LSS = Numericc::Solutions::LinearSystemSolution<double>;
    using PIS = Numericc::Solutions::PowerIterationSolution<double>;

    using namespace Numericc::RootFinding;

    RFS fixed = FixedPoint(g, -0.5, 100, 5);
    RFS bisection = Bisection(f1, 0, 2, 5, 5);
    RFS newton = NewtonRaphson(f1, df1, 1, 10, 5);
```
