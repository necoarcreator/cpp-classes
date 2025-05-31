#ifndef FILEREADCALC
#define FILEREADCALC

#include "CalculatorCommon.h"
#include <string>
#include <vector>
#include <map>
class FileReadCalc : public CalculatorCommon
{

protected:
	std::vector<std::vector<std::string>> originalString;
	std::map<const unsigned int, std::string> sequence;
	const std::string listOfCommands;
	long double result;
	unsigned int numRequests;
	unsigned int order;
	const int argc;
	std::vector<std::string> nameFiles;
	size_t numFiles;

	size_t movingOfIterator;

public:

	FileReadCalc(int _argc, std::vector<std::string> input, bool isWritingLog);
	~FileReadCalc();

	void const startCount() override;

	bool const check(char symb) override final;

	void const proceed(std::vector<std::string> input) override final;

	void const printSequence() override;

};

#endif