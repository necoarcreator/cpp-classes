#include <iostream>
#include <vector>
#include <string>
#include "Helper.h"
#include "StandartCalc.h"
#include "FileWriteCalc.h"
#include "FileReadWriteCalc.h"
#include "FileReadCalc.h"

using namespace std;

bool check(char symb)
{
	if ((static_cast<int>(symb) > 47) and (static_cast<int>(symb) < 58) or symb == '.' or symb == '(' or symb == ')')
	{
		return true;
	}
	return false;
}
vector<string> returnString(int argc, char** argv)
{
	vector<string> arggv;
	for (size_t i = 0; i < argc; i++)
	{
		string argI = argv[i];
		arggv.push_back(argI);
	}
	return arggv;
}
void checkCase(int _argc, std::vector<std::string> _input)
{
	size_t movingOfIterator = 0, movingIteratorForFile = 0;

	bool isWritingLog = false;

	if ((_input[1] == "-h") or (_input[1] == "-h "))
	{
		Helper A;
		return;
	}
	else if ((_input[1] == "calc") or (_input[1] == "calc ") or (_input[1] == " calc"))
	{

		if ((_input[2] == "-l") or (_input[2] == "-l ") or (_input[2] == " -l"))
		{
			isWritingLog = true;
			movingOfIterator = 1;
			movingIteratorForFile = 1;
			if (_argc > 5)
			{
				if ((_input[3] == "log" or _input[3] == " log" or _input[3] == "log ") and (not check(_input[5][0])))
				{

					vector<string> inputEff;

					inputEff.insert(inputEff.begin(), _input.begin(), _input.begin() + 3);
					inputEff.insert(inputEff.begin() + 3, _input.begin() + 5, _input.end());
					FileReadWriteCalc B(_argc, inputEff, isWritingLog, _input[4]);
					B.startCount();
					inputEff.clear();
					return;
				}
			}
			if (_argc > 4)
			{
				if (_input[3] == "log" or _input[3] == " log" or _input[3] == "log ")
				{

					vector<string> inputEff;

					inputEff.insert(inputEff.begin(), _input.begin(), _input.begin() + 3);
					inputEff.insert(inputEff.begin() + 3, _input.begin() + 5, _input.end());
					FileWriteCalc C(_argc, inputEff, isWritingLog, _input[4]);
					C.startCount();
					inputEff.clear();
					return;
				}
			}
		}
		if (_argc > 2)
		{
			if (not check(_input[2 + movingIteratorForFile][0]))
			{

				FileReadCalc D(_argc, _input, isWritingLog);
				D.startCount();
				return;
			}

		}


		StandartCalc E(_argc, _input, isWritingLog);
		E.startCount();
		return;
	}
	else
	{
		cerr << "Error! No command inserted. Try -h to open list of commands.\n";
		return;
	}
}

int main(int argc, char** argv)
{
	//const vector<string> ar{ "puk.exe", "calc", "exampleRequest.txt"};
	//int arg = ar.size();
	//cout << ar[1] << "   " << ar[2];
	checkCase(returnString(argc, argv).size(), returnString(argc, argv));
}