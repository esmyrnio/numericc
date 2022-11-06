/*******************************************************************
 *** numericc header file:
 ***
 *** ODE solver functions.
 ***
 *** Author: Evangelos Smyrniotis, November 2022
 ***
 *** The MIT License (MIT)
 *** Copyright (c) 2022 Evangelos Smyrniotis
 *** See LICENSE.txt
 *******************************************************************/
#pragma once

#include <iostream>

#include "../misc/Solutions.h"
#include "../misc/Massert.h"

namespace ODE
{
	using db = double;

	Solutions::OdeSolution Euler(db x0, db y0, db const step, db const x, db(*ODE)(db x, db y));

	namespace RK
	{
		db RK_one(db x, db y, db(*ODE)(db x, db y));
		db RK_two(db x, db y, db h, db(*ODE)(db x, db y));
		db RK_three(db x, db y, db h, db(*ODE)(db x, db y));
		db RK_four(db x, db y, db h, db(*ODE)(db x, db y));
		Solutions::OdeSolution RK4(db x0, db y0, db const step, db const x, db(*ODE)(db x, db y));
	}
}

#include "Differential.ipp"