/*******************************************************************
 *** numericc implementation file:
 ***
 *** Linear system solver and power-iteration method.
 ***
 *** Author: Evangelos Smyrniotis, November 2022
 ***
 *** The MIT License (MIT)
 *** Copyright (c) 2022 Evangelos Smyrniotis
 *** See LICENSE.txt
 *******************************************************************/
namespace LinearAlgebra
{
	/*
	  LinearSystemSolver solves a linear system of the form Ax = B,
	  via Gaussian elimination, e.g triangulation of the
	  augmented matrix A|B.

	  @ A_ : Coefficient matrix of the system A @ x = B.
	  @ B_ : Constant matrix of the system A @ x = B.

	  @ Returns a struct containing the solution matrix X
		as well as the A and B components of
		the A|B resulting triangular augmented matrix.
	*/

	template<typename T>
	Solutions::LinearSystemSolution<double> LinearSystemSolver(Super::Matrix<T> const& A_, Super::Matrix<T> const& B_)
	{
		Assert(A_.rows() == A_.cols(), "Linear system solver failed. Non-square coefficient matrix was given as input.");
		Assert(B_.rows() == 1 || B_.cols() == 1, "LinearSystemSolver failed. Non-vector constant matrix was given as input.");

		auto A = A_.Convert(); // to double type
		auto B = B_.rows() == 1 ? B_.Convert().Transpose() : B_.Convert(); // to double type

		auto N = A.rows();
		for (size_t k = 0; k < N - 1; ++k)
		{
			auto diag = std::abs(A[k][k]);
			if (!diag)
			{
				for (size_t i = k + 1; k < N; ++i)
				{

					bool pivot = std::abs(A[i][k]) > std::abs(A[k][k]);
					if (pivot)
					{
						std::vector<double> tempA;
						double tempB;

						for (size_t j = 0; j < N; j++)
						{
							tempA.push_back(A[i][j]);
							A[i][j] = A[k][j];
							A[k][j] = tempA.back();
						}
						tempB = B[k][0];
						B[k][0] = B[i][0];
						B[i][0] = tempB;
						break;
					}
				}
			}
			for (size_t i = k + 1; i < N; ++i)
			{
				if (!A[i][k]) continue;

				auto factor = A[k][k] / A[i][k];
				for (size_t j = k; j < N; j++)
				{
					A[i][j] = A[k][j] - A[i][j] * factor;
				}
				B[i][0] = B[k][0] - B[i][0] * factor;
			}
		}

		Super::Matrix<double> result(1, N);
		result[0][N - 1] = B[N - 1][0] / A[N - 1][N - 1];

		for (int i = N - 2; i >= 0; i--)
		{
			double s = 0;
			for (size_t j = i + 1; j < N; j++) s += A[i][j] * result[0][j];
			result[0][i] = (B[i][0] - s) / A[i][i];
		}

		Solutions::LinearSystemSolution<double> solution;

		solution.A = A;
		solution.B = B;
		solution.X = result.Transpose();

		return solution;
	}

	/*
	  PowerIteration computes the largest eigenvalue of a square matrix
	  as well as the corresponding eigenvector.

	  @ A_ : Coefficient matrix of the system A @ x = B.
	  @ B_ : Constant matrix of the system A @ x = B.

	  @ Returns a struct containing the largest eigenvalue of matrix A
		and the corresponding eigenvector.
	*/

	template<typename T>
	Solutions::PowerIterationSolution<double> PowerIteration(Super::Matrix<T> const& A_, size_t iterations)
	{
		Assert(A_.rows() == A_.cols(), "Power iteration method failed. Non-square coefficient matrix was given as input.");

		auto A = A_.Convert();
		std::default_random_engine generator;
		std::normal_distribution<double> distribution(0.0, 1.0);

		auto N = A.rows();

		std::vector<double> pseudoEigenVals;
		Super::Matrix<double> randMat(1, N);
		for (size_t i = 0; i < N; ++i) randMat[0][i] = distribution(generator);

		randMat = randMat.Transpose();

		size_t iter = 0;
		while (iter < iterations)
		{
			randMat = A * randMat;

			if (iter >= (iterations - 2))
				pseudoEigenVals.push_back(randMat.AbsoluteMax());
			iter++;
		}

		auto trandMat = randMat.Transpose();
		auto n_norm = trandMat * randMat;
		auto norm = sqrt(n_norm[0][0]);
		randMat = randMat / norm;

		Solutions::PowerIterationSolution<double> solution;

		solution.eigenvalue = pseudoEigenVals[1] / pseudoEigenVals[0];
		solution.eigenvector = randMat;

		return solution;
	}
}