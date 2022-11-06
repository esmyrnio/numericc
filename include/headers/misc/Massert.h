/*******************************************************************
 *** numericc header file:
 *** 
 *** Custom assert function and macro.
 *** 
 *** Author: Evangelos Smyrniotis, November 2022
 *** 
 *** The MIT License (MIT)
 *** Copyright (c) 2022 Evangelos Smyrniotis
 *** See LICENSE.txt
 *******************************************************************/
#pragma once

#include <iostream>

void _Assert(const char* condition_str, bool condition, const char* file, int line, const char* message);

#define Assert(condition, message) _Assert(#condition, condition, __FILE__, __LINE__, message)

#include "Massert.ipp"