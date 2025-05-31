#ifndef CALCCOMMON
#define CALCCOMMON

#include <string>
#include <vector>

class CalculatorCommon
{
public:

	virtual void const startCount() = 0;

	virtual bool const check(char symb) = 0;

	virtual void const proceed(std::vector<std::string> input) = 0;

	virtual void const printSequence() = 0;

};


#endif