#include "CrationalExeption.h"
#include <iostream>
#include <string>
using namespace std;

CrationalExeption::CrationalExeption(const char divZero[]) : whatToCout{ divZero } {};

CrationalExeption::CrationalExeption(float overFlow) : whatToCout { "Error! Your number is more than it's possible for float type!" } {};

const char* CrationalExeption::getError()
{
	return whatToCout;
}