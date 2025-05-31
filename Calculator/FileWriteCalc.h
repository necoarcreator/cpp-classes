#ifndef FILEWRITECALC
#define FILEWRITECALC

#include "StandartCalc.h"


class FileWriteCalc : public StandartCalc
{
protected:
	std::string nameFile;
public:
	void const printSequence() override final;
	void const startCount() override final;
	FileWriteCalc(int _argc, std::vector<std::string> input, bool isWritingLog, std::string _nameFile);
	~FileWriteCalc();
};

#endif