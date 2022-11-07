/*******************************************************************
 *** numericc implementation file:
 ***
 *** Root finding routines.
 ***
 *** Author: Evangelos Smyrniotis, November 2022
 ***
 *** The MIT License (MIT)
 *** Copyright (c) 2022 Evangelos Smyrniotis
 *** See LICENSE.txt
 *******************************************************************/
namespace RootFinding
{
    /*
      Bisection method implementation.

      @ f : function for which the root is to be found.
      @ a : left limit of root bracketing interval.
      @ b : right limit of root bracketing interval.
      @ maxIter : maximum number of iterations allowed.
      @ scarboroughCrit: Scarborough criterion's parameter for tracking convergence of the solution.

      @ Returns a struct containing the root solution, the total number of iterations,
        and a boolean convergence check.
    */

    Solutions::RootFindingSolution Bisection(db(*f)(db x), db a, db b, size_t maxIter, size_t scarboroughCrit)
    {
        Assert(f(a) * f(b) < 0, "BisectionMethod: requires root bracketing.");
        Solutions::RootFindingSolution solution(0, 0, false);
        db currmid, prevmid, fmid, checkLeft, checkRight;
        for (int n = 0; n < maxIter; n++)
        {
            currmid = (a + b) / 2;
            fmid = f(currmid);
            checkLeft = f(a) * fmid;
            checkRight = f(b) * fmid;
            if (checkLeft < 0) b = currmid;
            else if (checkRight < 0) a = currmid;
            else if (fmid == 0)
            {
                solution(currmid, n, true);
                return solution;
            }
            else
            {
                std::cerr << "Bisection method failed.\n";
                exit(0);
            }

            if (n >= 2 && std::abs((currmid - prevmid) / currmid) < 0.5 * pow(10, 2 - static_cast<db>(scarboroughCrit)))
            {
                solution(currmid, n, true);
                return solution;
            }
            prevmid = currmid;
        }
        //Assert(solution.converged, "Bisection method failed to converge for the given input parameters.");
        return solution;
    }

    /*
      Newton-Raphson method implementation.

      @ f : function for which the root is to be found.
      @ startingPoint : starting guess solution for the root.
      @ maxIter : maximum number of iterations allowed.
      @ scarboroughCrit: Scarborough criterion's parameter for tracking convergence of the solution.

      @ Returns a struct containing the root solution, the total number of iterations,
        and a boolean convergence check.
    */

    Solutions::RootFindingSolution NewtonRaphson(db(*f)(db x), db(*df)(db x), const db startingPoint, size_t maxIter, size_t scarboroughCrit)
    {
        Assert(df(startingPoint) != 0, "NewtonRaphson: requires non-zero first order derivative f'(x) at x = guess.");

        size_t n = 0;
        auto r = f(startingPoint) / df(startingPoint);
        auto check = 0.5 * pow(10, 2 - static_cast<db>(scarboroughCrit));
        auto currentPoint = startingPoint;
        while (std::abs(r) >= check && n <= maxIter)
        {
            currentPoint -= r;
            n++;
            Assert(df(currentPoint) != 0, "Division by zero detected. Convergence failed.");
            r = f(currentPoint) / df(currentPoint);
        }
        Assert(n < maxIter, "Newton-Raphson method failed to converge for the given input parameters.");
        Solutions::RootFindingSolution solution = { currentPoint, n , true };
        return solution;
    }

    /*
      Fixed-point method implementation.

      @ f : g(x) function (g(x)=x) for initial function for which the root is to be found.
      @ startingPoint : starting guess solution for the root.
      @ maxIter : maximum number of iterations allowed.
      @ scarboroughCrit: Scarborough criterion's parameter for tracking convergence of the solution.

      @ Returns a struct containing the root solution and the total number of iterations.
    */

    Solutions::RootFindingSolution FixedPoint(db(*f)(db x), const db startingPoint, size_t maxIter, size_t scarboroughCrit)
    {
        db prev, curr;
        prev = startingPoint;
        size_t n = 1;
        Solutions::RootFindingSolution solution(0, 0, false);
        if (!f(prev)) return solution(prev, 1, true);
        auto check = 0.5 * pow(10, 2 - static_cast<db>(scarboroughCrit));
        while (n < maxIter)
        {
            curr = f(prev);

            if (curr != 0 && std::abs((curr - prev) / curr) < check) return solution(curr, n, true);
            n++;
            prev = curr;
        }
        //Assert(n < maxIter, "Fixed-point iteration failed to converge for the given input parameters.");
        return solution;
    }
}
