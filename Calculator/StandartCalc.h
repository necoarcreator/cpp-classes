#ifndef STANDARTCALC
#define STANDARTCALC

#include "CalculatorCommon.h"
#include <string>
#include <vector>
#include <map>

class StandartCalc : public CalculatorCommon
{
protected:
	std::vector<std::vector<std::string>> originalString;
	std::map<const unsigned int, std::string> sequence;
	const std::string listOfCommands;
	long double result;
	unsigned int numRequests;
	unsigned int order;
	const int argc;

	size_t movingOfIterator;

public:
	StandartCalc();
	StandartCalc(int _argc, std::vector<std::string> input, bool isWritingLog);
	~StandartCalc();

	void const startCount() override;

	bool const check(char symb) override final;

	void const proceed(std::vector<std::string> input) override final;

	void const printSequence() override;

};

#endif