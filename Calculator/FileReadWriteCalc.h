#ifndef FILEREADWRITE
#define FILEREADWRITE

#include "FileReadCalc.h"
#include <fstream>
#include <iostream>

#include <vector>
#include <string>
using namespace std;

class FileReadWriteCalc : public FileReadCalc
{
protected:

	std::string nameFile;

public:

	void const printSequence() override final;
	void const startCount() override final;
	FileReadWriteCalc(int _argc, std::vector<std::string> input, bool isWritingLog, std::string _nameFile);
	~FileReadWriteCalc();
};

#endif