#include "FileWriteCalc.h"
#include <fstream>
#include <iostream>

#include <vector>
#include <string>

using namespace std;

FileWriteCalc::FileWriteCalc(int _argc, vector<string> input, bool isWritingLog, string _nameFile) : StandartCalc(_argc - 2, input, isWritingLog), nameFile(_nameFile)
{}

void const FileWriteCalc::printSequence()
{
	ofstream out;
	out.open(nameFile, fstream::app);
	if (out.is_open())
	{
		out << "********************************" << endl;
		for (auto it = sequence.begin(); it != sequence.end(); ++it)
		{
			out << "operation number " << it->first << ": " << it->second << endl;
			out << "********************************" << endl;
		}
	}
	else
	{
		cerr << "Error! Can't open a file for output!\n";
		exit(1);
	}
	out.close();

}
void const FileWriteCalc::startCount()
{
	ofstream out0(nameFile);
	out0.close();

	for (size_t i = 0; i < numRequests; i++)
	{
		proceed(originalString[i]);
		ofstream out;
		out.open(nameFile, fstream::app);
		if (out.is_open())
		{
			out << result << endl;
			out << endl;
			result = 0;
			result = 0;
		}
		else
		{
			cerr << "Error! Can't open a file for output!\n";
			exit(1);
		}
		out.close();

		if (movingOfIterator == 1)
		{
			printSequence();
			sequence.clear();
			order = 1;
		}
	}
}
FileWriteCalc::~FileWriteCalc()
{

	cout << "Please don't look here, the result is counted in a file.\n";
	originalString.clear();
	sequence.clear();
	nameFile.clear();
}