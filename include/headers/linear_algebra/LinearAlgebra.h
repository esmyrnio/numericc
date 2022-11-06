/*******************************************************************
 *** numericc header file:
 ***
 *** Linear system solver and power-iteration method.
 ***
 *** Author: Evangelos Smyrniotis, November 2022
 ***
 *** The MIT License (MIT)
 *** Copyright (c) 2022 Evangelos Smyrniotis
 *** See LICENSE.txt
 *******************************************************************/
#pragma once

#include <random>
#include <cmath>
#include "Matrix.h"

namespace LinearAlgebra
{
	template<typename T>
	Solutions::LinearSystemSolution<double> LinearSystemSolver(Super::Matrix<T> const& A_, Super::Matrix<T> const& B_);

	template<typename T>
	Solutions::PowerIterationSolution<double> PowerIteration(Super::Matrix<T> const& A_, size_t iterations);
}

#include "LinearAlgebra.ipp"
