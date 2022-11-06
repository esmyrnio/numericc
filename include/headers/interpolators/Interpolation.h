/*******************************************************************
 *** numericc header file:
 ***
 *** Interpolation routines.
 ***
 *** Author: Evangelos Smyrniotis, November 2022
 ***
 *** The MIT License (MIT)
 *** Copyright (c) 2022 Evangelos Smyrniotis
 *** See LICENSE.txt
 *******************************************************************/
#pragma once

#include <vector>
#include <cmath>

#include "../misc/Massert.h"

namespace Interpolators
{
	using db = double;
	db Lagrange(const db x, const std::vector<db>& XPoints, const std::vector<db>& YPoints);

	namespace Newton
	{
		size_t factorial(size_t x);
		db ufunc(db u, size_t n);
		db NewtonGregoryFD(const db x, std::vector<db> const& XPoints, std::vector<db> const& YPoints);
		db NewtonGregoryBD(const db x, std::vector<db> const& XPoints, std::vector<db> const& YPoints);
	}
}

#include "Interpolation.ipp"