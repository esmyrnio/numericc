/*******************************************************************
 *** numericc header file:
 ***
 *** Root finding routines.
 ***
 *** Author: Evangelos Smyrniotis, November 2022
 ***
 *** The MIT License (MIT)
 *** Copyright (c) 2022 Evangelos Smyrniotis
 *** See LICENSE.txt
 *******************************************************************/
#pragma once

#include "../misc/Solutions.h"
#include "../misc/Massert.h"

namespace RootFinding
{
    using db = double;

    Solutions::RootFindingSolution Bisection(db(*f)(db x), db a, db b, size_t maxIter, size_t scarboroughCrit);
    Solutions::RootFindingSolution NewtonRaphson(db(*f)(db x), db(*df)(db x), const db startingPoint, size_t maxIter, size_t scarboroughCrit);
    Solutions::RootFindingSolution FixedPoint(db(*f)(db x), const db startingPoint, size_t maxIter, size_t scarboroughCrit);
}

#include "RootFinding.ipp"