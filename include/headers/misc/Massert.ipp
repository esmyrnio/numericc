/*******************************************************************
 *** numericc implementation file:
 ***
 *** Custom assert function.
 ***
 *** Author: Evangelos Smyrniotis, November 2022
 ***
 *** The MIT License (MIT)
 *** Copyright (c) 2022 Evangelos Smyrniotis
 *** See LICENSE.txt
 *******************************************************************/
void _Assert(const char* condition_str, bool condition, const char* file, int line, const char* message)
{
    if (!condition)
    {
        std::cerr << "Assert failed:\t" << message << "\n"
            << "Expected:\t" << condition_str << "\n"
            << "Source:\t\t" << file << ", line " << line << "\n";
        abort();
    }
}