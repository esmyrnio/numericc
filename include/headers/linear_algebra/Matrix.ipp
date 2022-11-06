/*******************************************************************
 *** numericc implementation file:
 ***
 *** Matrix class.
 ***
 *** Author: Evangelos Smyrniotis, November 2022
 ***
 *** The MIT License (MIT)
 *** Copyright (c) 2022 Evangelos Smyrniotis
 *** See LICENSE.txt
 *******************************************************************/
namespace Super
{

	// **************************************************************************************************************
	//	@ Constructors
	// **************************************************************************************************************

	template<typename T>
	Matrix<T>::Matrix() { _inner.emplace_back(Row(1)); }

	template<typename T>
	Matrix<T>::Matrix(size_t N) { _inner.emplace_back(Row(N)); }

	template<typename T>
	Matrix<T>::Matrix(size_t r, size_t c) : _inner(r, Row(c, static_cast<T>(0))) {}

	template<typename T>
	Matrix<T>::Matrix(Init const& Input) { _inner.emplace_back(Input); }

	template<typename T>
	Matrix<T>::Matrix(Init2d const& Input)
	{
		for (auto&& r : Input) _inner.emplace_back(r);
	}

	template<typename T>
	Matrix<T>::Matrix(std::vector<Row> const& Input)
	{
		_inner.reserve(Input.size());
		for (auto&& r : Input) _inner.emplace_back(r);
	}

	template<typename T>
	Matrix<T>::Matrix(Row const& Input)
	{
		_inner.emplace_back(Input);
	}

	// **************************************************************************************************************
	//	@ Operators
	// **************************************************************************************************************

	template<typename T>
	inline typename Matrix<T>::Row& Matrix<T>::operator[](size_t a)
	{
		return _inner[a];
	}

	template<typename T>
	inline typename const Matrix<T>::Row const& Matrix<T>::operator[](size_t a) const
	{
		return _inner[a];
	}

	template<typename T>
	Matrix<T> operator+(Matrix<T> const& A, Matrix<T> const& B)
	{
		Assert(A.cols == B.cols && A.rows() == B.rows(), "Matrix addition failed. Matrices must have the same number of rows and columns.");
		std::vector<std::vector<T>> starting(A.rows(), std::vector<T>(A.cols, static_cast<T>(0)));
		Matrix<T> result(starting);

		for (size_t i = 0; i < A.rows(); ++i)
			for (size_t j = 0; j < A.cols(); ++j)
			{
				result[i][j] = A[i][j] + B[i][j];
			}
		return result;
	}

	template<typename T>
	Matrix<T> operator-(Matrix<T> const& A, Matrix<T> const& B)
	{
		Assert(A.cols() == B.cols() && A.rows() == B.rows(), "Matrix substraction failed. Matrices must have the same number of rows and columns.");
		std::vector<std::vector<T>> starting(A.rows(), std::vector<T>(A.cols(), static_cast<T>(0)));
		Matrix<T> result(starting);

		for (size_t i = 0; i < A.rows(); ++i)
			for (size_t j = 0; j < A.cols(); ++j)
			{
				result[i][j] = A[i][j] - B[i][j];
			}
		return result;
	}

	template<typename T>
	Matrix<T> operator*(Matrix<T> const& A, Matrix<T> const& B)
	{

		Assert(A.cols() == B.rows(), "Matrix multiplication failed. Column number of matrix A not equal to row number of matrix B");
		std::vector<std::vector<T>> starting(A.rows(), std::vector<T>(B.cols(), static_cast<T>(0)));
		Matrix<T> result(starting);
		for (size_t i = 0; i < A.rows(); ++i)
			for (size_t j = 0; j < B.cols(); ++j)
			{
				result[i][j] = 0;
				for (size_t k = 0; k < B.rows(); ++k) result[i][j] += A[i][k] * B[k][j];
			}
		return result;
	}

	template<typename T>
	Matrix<T> operator*(Matrix<T> const& A, T const& val)
	{

		Matrix<T> result = A;
		for (size_t i = 0; i < A.rows(); ++i)
			for (size_t j = 0; j < A.cols(); ++j)
			{
				result[i][j] *= val;
			}
		return result;
	}

	template<typename T>
	Matrix<T> operator/(Matrix<T> const& A, T const& val)
	{

		Matrix<T> result = A;
		for (size_t i = 0; i < A.rows(); ++i)
			for (size_t j = 0; j < A.cols(); ++j)
			{
				result[i][j] /= val;
			}
		return result;
	}
	template<typename T>
	bool operator==(Matrix<T> const& A, Matrix<T> const& B)
	{
		return A._inner == B._inner;
	}
	template<typename T>
	bool operator!=(Matrix<T> const& A, Matrix<T> const& B)
	{
		return !(A._inner == B._inner);
	}
	template<typename T>
	std::ostream& operator<<(std::ostream& os, Matrix<T> const& m)
	{
		os << "\n[";
		for (size_t i = 0; i < m.rows(); i++)
		{
			for (size_t j = 0; j < m.cols(); j++)
			{
				auto val = (j < m.cols() - 1 ? "\t" : " ");
				os << (j == 0 && i != 0 && m[i][j] >= 0 ? "  " : " ") << m[i][j] << val;
			}
			auto val = (i < m.rows() - 1 ? "\n" : "]\n");
			os << val;
		}
		return os;
	}

	// **************************************************************************************************************
	//	@ Routines
	// **************************************************************************************************************

	template<typename T>
	Solutions::MinMaxSolution<T> Matrix<T>::MinMax() const
	{
		auto self = *this;

		auto Nx = self.cols();
		auto Ny = self.rows();

		T max = self[0][0];
		T min = self[0][0];

		for (size_t i = 0; i < Ny; ++i)
			for (size_t j = 0; j < Nx; ++j)
			{
				if (self[i][j] > max) max = self[i][j];
				if (self[i][j] < min) min = self[i][j];
			}

		Solutions::MinMaxSolution<T> minmax;
		minmax.min = min;
		minmax.max = max;
		return minmax;
	}

	template<typename T>
	T Matrix<T>::AbsoluteMax() const
	{
		auto self = *this;

		auto Nx = self.cols();
		auto Ny = self.rows();

		T max = self[0][0];

		for (size_t i = 0; i < Ny; ++i)
			for (size_t j = 0; j < Nx; ++j)
			{
				if (std::abs(self[i][j]) > max) max = std::abs(self[i][j]);
			}

		return max;
	}

	template<typename T>
	Matrix<T> Matrix<T>::Transpose() const
	{
		auto self = *this;
		Matrix<T> result(self.cols(), self.rows());
		for (size_t i = 0; i < result.rows(); i++)
			for (size_t j = 0; j < result.cols(); j++)
				result[i][j] = self[j][i];

		return result;
	}

	template<typename T>
	Matrix<T> Matrix<T>::Cofactor(size_t line, size_t column, size_t size) const
	{
		auto self = *this;
		Assert(self.rows() == self.cols(), "Cofactor failed: non-square matrix.");

		Matrix<T> result = self;
		size_t N = self.m.rows();
		size_t i = 0, j = 0;
		for (size_t r = 0; r < size; ++r) {
			for (size_t c = 0; c < size; ++c) {
				if (r != line && c != column) {
					result[i][j++] = self[r][c];
					if (j == size - 1) {
						j = 0;
						i++;
					}
				}
			}
		}
		return result;
	}
	template<typename T>
	T Matrix<T>::Determinant()
	{
		return Determinant(this->m.rows());
	}
	template<typename T>
	T Matrix<T>::Determinant(size_t size) const
	{
		auto self = *this;
		Assert(self.rows() == self.cols(), "Determinant routine failed: non-square matrix.");
		T det = 0;
		if (!(size - 1)) return self[0][0];
		Mat cofactors;
		int sign = 1;
		for (size_t c = 0; c < size; c++) {
			cofactors = self.Cofactor(0, c, size);
			det += sign * self[0][c] * cofactors.Determinant(size - 1);
			sign = -sign;
		}
		return det;
	}

	template<typename T>
	Matrix<T> Matrix<T>::Adjoint() const
	{
		auto self = *this;
		Assert(self.rows() == self.cols(), "Adjoint routine failed: non-square matrix.");
		auto N = self.m.rows();

		Mat result = self;

		if (N == 1)
		{
			result[0][0] = 1;
			return result;
		}

		int sign = 1;
		Mat cofactors;

		for (size_t i = 0; i < N; i++) {
			for (size_t j = 0; j < N; j++) {
				cofactors = self.Cofactor(i, j, N);
				sign = ((i + j) % 2 == 0) ? 1 : -1;
				result[j][i] = sign * cofactors.Determinant(N - 1);
			}
		}
		return result;
	}
	template<typename T>
	Matrix<double> Matrix<T>::Inverse() const
	{
		auto self = *this;
		Assert(self.rows() == self.cols(), "Inverse routine failed: non-square matrix.");

		auto N = self.m.rows();
		auto det = self.Determinant(N);
		Assert(det != 0, "Inverse routine failed. Determinant of matrix is zero.");
		Matrix<double> result(self.rows(), self.cols());
		Mat adjoint = self.Adjoint();
		for (size_t i = 0; i < N; ++i)
			for (size_t j = 0; j < N; ++j)
				result[i][j] = static_cast<double>(adjoint[i][j]) / det;
		return result;
	}

	template<typename T>
	Matrix<double> Matrix<T>::Convert() const
	{
		auto self = *this;
		Matrix<double> db(self.rows(), self.cols());
		for (size_t i = 0; i < self.rows(); ++i)
		{
			for (size_t j = 0; j < self.cols(); ++j)
			{
				db[i][j] = (double)self[i][j];
			}
		}
		return db;
	}
}