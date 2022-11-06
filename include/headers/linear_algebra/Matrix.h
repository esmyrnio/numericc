/*******************************************************************
 *** numericc header file:
 ***
 *** Matrix class.
 ***
 *** Author: Evangelos Smyrniotis, November 2022
 ***
 *** The MIT License (MIT)
 *** Copyright (c) 2022 Evangelos Smyrniotis
 *** See LICENSE.txt
 *******************************************************************/
#pragma once

#include <vector>

#include "../misc/Solutions.h"
#include "../misc/Massert.h"

namespace Super
{
	template<typename T>
	class Matrix
	{
		private:
		
			using Row = typename std::vector<T>;
			using Init = typename std::initializer_list<T>;
			using Init2d = typename std::initializer_list<std::initializer_list<T>>;
			using Mat = typename Matrix<T>;

			std::vector<Row> _inner;

		public:

			// *******************************************************
			//	@ Constructors
			// *******************************************************

			Matrix();
			Matrix(size_t N);
			Matrix(size_t r, size_t c);
			Matrix(Init const& Input);
			Matrix(Init2d const& Input);
			Matrix(std::vector<Row> const& Input);
			Matrix(Row const& Input);

			// *******************************************************
			//	@ Dimensions
			// *******************************************************

			constexpr size_t rows() const { return _inner.size(); };
			constexpr size_t cols() const { return _inner[0].size(); };
			constexpr size_t m_size() const { return rows() * cols(); };

			// *******************************************************
			//	@ Routines
			// *******************************************************

			T Determinant();
			T Determinant(size_t size) const;
			Mat Transpose() const;
			Mat Cofactor(size_t line, size_t column, size_t size) const;
			Mat Adjoint() const;
			Matrix<double> Inverse() const;

			Solutions::MinMaxSolution<T> MinMax() const;
			T AbsoluteMax() const;

			Matrix<double> Convert() const;
		
			// *******************************************************
			//	@ Operator overloading (internal/friend)
			// *******************************************************

			Row& operator[](size_t a);
			Row const& operator[](size_t a) const;
			template<typename U>
			friend std::ostream& operator<<(std::ostream& os, Matrix<U> const& m);
			template<typename U>
			friend bool operator==(Matrix<U> const& A, Matrix<U> const& B);
			template<typename U>
			friend bool operator!=(Matrix<U> const& A, Matrix<U> const& B);
	};

	// *******************************************************
	//	@ Operator overloading (external)
	// *******************************************************

	template<typename T>
	Matrix<T> operator+(Matrix<T> const& A, Matrix<T> const& B);
	template<typename T>
	Matrix<T> operator-(Matrix<T> const& A, Matrix<T> const& B);
	template<typename T>
	Matrix<T> operator*(Matrix<T> const& A, Matrix<T> const& B);
	template<typename T>
	Matrix<T> operator*(Matrix<T> const& A, T const& val);
	template<typename T>
	Matrix<T> operator/(Matrix<T> const& A, T const& val);
}

#include "Matrix.ipp"
