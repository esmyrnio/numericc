/*******************************************************************
 *** numericc implementation file:
 ***
 *** Interpolation routines.
 ***
 *** Author: Evangelos Smyrniotis, November 2022
 ***
 *** The MIT License (MIT)
 *** Copyright (c) 2022 Evangelos Smyrniotis
 *** See LICENSE.txt
 *******************************************************************/
namespace Interpolators
{
	/*
	  Lagrange interpolation method implementation.

	  @ x : point to evaluate the intepolation.
	  @ XPoints : vector of independent variable points.
	  @ YPoints : vector of dependent variable points.

	  @ Returns the interpolated value of y at given x.
	*/

	db Lagrange(const db x, const std::vector<db>& XPoints, const std::vector<db>& YPoints)
	{
		size_t n = XPoints.size();
		db s = 0;
		for (size_t i = 0; i < n; i++)
		{
			db pol = 1;
			for (size_t j = 0; j < n; j++)
			{
				if (j != i) pol *= (x - XPoints[j]) / (XPoints[i] - XPoints[j]);
			}
			s += YPoints[i] * pol;
		}
		return s;
	}

	namespace Newton
	{
		// @Returns the factorial of a given number.
		size_t factorial(size_t x)
		{
			size_t temp = 1;
			for (size_t i = 2; i <= x; i++) temp *= i;
			return temp;
		}

		// @Returns the u functions in Newton-Gregory interpolation formula.
		db ufunc(db u, size_t n)
		{
			auto x = u;
			for (size_t i = 1; i < n; i++) x *= u - 1;
			return x;
		}

		/*
		  Newton-Gregory forward interpolation method implementation.

		  @ x : point to evaluate the intepolation.
		  @ XPoints : vector of independent variable points.
		  @ YPoints : vector of dependent variable points.

		  @ Returns the interpolated value of y at given x.
		*/

		db NewtonGregoryFD(const db x, std::vector<db> const& XPoints, std::vector<db> const& YPoints)
		{
			auto Xsize = XPoints.size();
			auto Ysize = YPoints.size();
			auto N = Xsize;

			Assert(Xsize == Ysize, "NewtonGregory: requires equal number of X and Y input values.");

			std::vector<std::vector<db>> y(N, std::vector<db>(N));

			for (size_t k = 0; k < N; k++) y[k][0] = YPoints[k];

			for (size_t i = 1; i < N; i++) {
				for (size_t j = 0; j < N - i; j++)
					y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
			}

			auto s = y[0][0];
			auto u = (x - XPoints[0]) / (XPoints[1] - XPoints[0]);

			for (size_t i = 1; i < N; i++)
				s += (ufunc(u, i) * y[0][i]) / factorial(i);

			return s;
		};

		/*
		  Newton-Gregory backwards interpolation method implementation.

		  @ x : point to evaluate the intepolation.
		  @ XPoints : vector of independent variable points.
		  @ YPoints : vector of dependent variable points.

		  @ Returns the interpolated value of y at given x.
		*/

		db NewtonGregoryBD(const db x, std::vector<db> const& XPoints, std::vector<db> const& YPoints)
		{
			auto Xsize = XPoints.size();
			auto Ysize = YPoints.size();
			auto N = Xsize;

			Assert(Xsize == Ysize, "NewtonGregory: requires equal number of X and Y input values.");

			std::vector<std::vector<db>> y(N, std::vector<db>(N));

			for (size_t k = 0; k < N; k++) y[k][0] = YPoints[k];

			for (size_t i = 1; i < N; i++) {
				for (size_t j = N - 1; j >= 1; j--)
					y[j][i] = y[j][i - 1] - y[j - 1][i - 1];
			}

			auto s = y[N - 1][0];
			auto u = (x - XPoints[N - 1]) / (XPoints[1] - XPoints[0]);
			for (size_t i = 1; i < N; i++)
				s += (ufunc(u, i) * y[N - 1][i]) / factorial(i);

			return s;
		};
	}
}