/*******************************************************************
 *** numericc implementation file:
 ***
 *** ODE solver functions.
 ***
 *** Author: Evangelos Smyrniotis, November 2022
 ***
 *** The MIT License (MIT)
 *** Copyright (c) 2022 Evangelos Smyrniotis
 *** See LICENSE.txt
 *******************************************************************/
namespace ODE
{
	/*
	  Euler method implementation for solving first order ODEs.

	  @ x0 : independent variable initial value.
	  @ y0 : y(x0) initial condition.
	  @ step : step size.
	  @ x : independent variable.
	  @ ODE : right-hand side of ODE.

	  @ Returns a struct containing the value of computed y at given x,
	  as well as the vector solutions for x and y.
	*/

	Solutions::OdeSolution Euler(db x0, db y0, db const step, db const x, db(*ODE)(db x, db y))
	{
		Assert(step != 0, "Euler method failed: zero step size was given.");

		auto N = (x - x0) / step;

		std::vector<db> xs = { x0 };
		std::vector<db> ys = { y0 };
		db df, ynew;
		size_t k = 0;

		while (k < N)
		{
			df = ODE(x0, y0);
			ynew = y0 + step * df;
			y0 = ynew;
			x0 += step;
			xs.push_back(x0);
			ys.push_back(y0);
			k++;
		}

		Solutions::OdeSolution solution;
		solution.Xs = xs;
		solution.Ys = ys;
		solution.X = x;
		solution.Y = ys.back();
		return solution;
	}

	namespace RK
	{
		db RK_one(db x, db y, db(*ODE)(db x, db y)) { return ODE(x, y); }
		db RK_two(db x, db y, db h, db(*ODE)(db x, db y)) { return ODE(x + 0.5 * h, y + 0.5 * RK_one(x, y, ODE) * h); }
		db RK_three(db x, db y, db h, db(*ODE)(db x, db y)) { return ODE(x + 0.5 * h, y + 0.5 * RK_two(x, y, h, ODE) * h); }
		db RK_four(db x, db y, db h, db(*ODE)(db x, db y)) { return ODE(x + h, y + RK_three(x, y, h, ODE) * h); }

		/*
		  Runge-Kutta 4th-order method implementation for solving 1st order ODEs.

		  @ x0 : independent variable initial value.
		  @ y0 : y(x0) initial condition.
		  @ step : step size.
		  @ x : independent variable.
		  @ ODE : right-hand side of ODE.

		  @ Returns a struct containing the value of computed y at given x,
		  as well as the vector solutions for x and y.
		*/

		Solutions::OdeSolution RK4(db x0, db y0, db const step, db const x, db(*ODE)(db x, db y))
		{
			Assert(step != 0, "Euler method failed: zero step size was given.");
			auto N = (x - x0) / step;

			std::vector<db> xs = { x0 };
			std::vector<db> ys = { y0 };
			db ynew;
			for (size_t i = 0; i < N; ++i)
			{
				ynew = y0 + (1.0 / 6.0) * (RK_one(xs[i], ys[i], ODE)
					+ 2.0 * RK_two(xs[i], ys[i], step, ODE) + 2.0 * RK_three(xs[i], ys[i], step, ODE)
					+ RK_four(xs[i], ys[i], step, ODE)) * step;
				y0 = ynew;
				x0 += step;
				xs.push_back(x0);
				ys.push_back(y0);
			}

			Solutions::OdeSolution solution;
			solution.Xs = xs;
			solution.Ys = ys;
			solution.X = x;
			solution.Y = ys.back();
			return solution;
		}
	}
}
