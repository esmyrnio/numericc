/*******************************************************************
 *** numericc header file:
 ***
 *** Solution structs for library routines.
 ***
 *** Author: Evangelos Smyrniotis, November 2022
 ***
 *** The MIT License (MIT)
 *** Copyright (c) 2022 Evangelos Smyrniotis
 *** See LICENSE.txt
 *******************************************************************/
#pragma once

#include <vector>

namespace Super
{
	template<typename T>
	class Matrix;
}

namespace Solutions
{
	struct RootFindingSolution
	{
		RootFindingSolution(double root, size_t iterations, bool converged) : root(root), iterations(iterations), converged(converged) {};
		RootFindingSolution operator () (double root, size_t iterations, bool converged)
		{
			RootFindingSolution solution = { root, iterations, converged };
			return solution;
		}
		double root;
		size_t iterations;
		bool converged;
	};

	template<typename T>
	struct LinearSystemSolution
	{
		Super::Matrix<T> X, A, B;
	};

	template<typename T>
	struct PowerIterationSolution
	{
		T eigenvalue;
		Super::Matrix<T> eigenvector;
	};

	template<typename T>
	struct MinMaxSolution
	{
		T min, max;
	};

	struct OdeSolution
	{
		std::vector<double> Xs, Ys;
		double X,Y;
	};
}