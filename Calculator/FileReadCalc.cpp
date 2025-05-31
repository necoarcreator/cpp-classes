#include "FileReadCalc.h"

#include <fstream>
#include <iostream>

#include <vector>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;

FileReadCalc::FileReadCalc(int _argc, vector<string> input, bool isWritingLog) : listOfCommands("+-*/%^"), argc(_argc), order(1), numRequests(0), numFiles(0)
{
	if (isWritingLog)
	{
		movingOfIterator = 1;
	}
	else
	{
		movingOfIterator = 0;
	}

	for (size_t i = 2 + movingOfIterator; i < _argc; i++)
	{
		nameFiles.push_back(input[i]);
		numFiles++;
	}

	for (size_t i = 0; i < numFiles; i++)
	{
		ifstream file(nameFiles[i], ios_base::in);
		string line;
		if (file.is_open())
		{
			char const delim = '\n';
			char const delim2 = '\"';
			while (getline(file, line))
			{
				string toParce;
				stringstream streamToParce(line);

				getline(streamToParce, toParce, delim2);
				getline(streamToParce, toParce, delim2);
				if (toParce == "")
				{
					file.close();
					break;
				}
				numRequests++;
				originalString.push_back({});
				originalString[numRequests - 1].push_back("");
				size_t p = 0; //счетчик длины числа
				size_t j = 0; //счетчик числа элементов в строке ввода
				for (size_t k = 0; k < toParce.size(); k++)
				{//

					if (check(toParce[k]))
					{
						p++;

					}

					if ((toParce[k] == ' ') and (check(toParce[k - 1])))
					{
						string sub = toParce.substr(k - p, p);

						originalString[numRequests - 1][j] += sub;
						originalString[numRequests - 1].push_back("");
						j++;
						p = 0;
					}
					else if (not check(toParce[k]))
					{

						for (auto x : listOfCommands)
						{
							if (x == toParce[k])
							{
								if (check(toParce[k - 1]))
								{
									string sub = toParce.substr(k - p, p);

									originalString[numRequests - 1][j] += sub;
									originalString[numRequests - 1].push_back("");
									j++;
									p = 0;
								}

								originalString[numRequests - 1][j] += toParce[k]; 
								originalString[numRequests - 1].push_back("");
								j++;

								if (check(toParce[k + 1]))
								{
									p++;
								}
								p = 0;
							}
						}

					}

					if (k == toParce.size() - 1)
					{
						string sub = toParce.substr(k - p + 1, p);

						originalString[numRequests - 1][j] += sub;

					}

				}
			}
		}
		else
		{
			cerr << "Error while opening the file to read!\n";
			exit(1);
		}
		file.close();
	}

}

void const FileReadCalc::printSequence()
{
	cout << "********************************" << endl;
	for (auto it = sequence.begin(); it != sequence.end(); ++it)
	{
		cout << "operation number " << it->first << ": " << it->second << endl;
		cout << "********************************" << endl;
		cout << endl;
	}
}

FileReadCalc::~FileReadCalc()
{
	originalString.clear();
	sequence.clear();
	nameFiles.clear();
}

void const FileReadCalc::startCount()
{
	for (size_t i = 0; i < numRequests; i++)
	{
		proceed(originalString[i]);
		cout << result << endl;
		result = 0;
		if (movingOfIterator == 1)
		{
			printSequence();
			sequence.clear();
			order = 1;
		}
	}
}
void const FileReadCalc::proceed(vector<string> input)
{

	size_t startBrackets;
	size_t endBrackets;
	for (size_t i = 0; i < input.size(); i++)
	{
		if (input[i][0] == '(')
		{
			startBrackets = i;
		}
		if (input[i][input[i].size() - 1] == ')')
		{
			endBrackets = i;
			vector<string> partOfInput;
			string inputWithoutStartBracket = input[startBrackets].substr(1, input[startBrackets].size() - 1);
			string inputWithoutEndBracket = input[endBrackets].substr(0, input[endBrackets].size() - 1);
			partOfInput.push_back(inputWithoutStartBracket);
			partOfInput.insert(partOfInput.begin() + 1, input.begin() + startBrackets + 1, input.begin() + endBrackets);
			partOfInput.push_back(inputWithoutEndBracket);

			proceed(partOfInput);
			input.erase(input.begin() + startBrackets, input.begin() + endBrackets + 1);
			input.insert(input.begin() + startBrackets, to_string(result));
		}


	}

	while (true)
	{

		string whatOperation;
		for (size_t i = 0; i < input.size(); i++)
		{
			long double localResult = 0;

			if ((input[i][0] == '*') or (input[i][0] == '/') or (input[i][0] == '%') or (input[i][0] == '^')) //сначала - деление, умножение, деление нацело, возведение в степень
			{
				long double a = 0.0;
				long double b = 0.0;
				size_t j = i - 1, k = i + 1;
				switch (input[i][0])
				{
				case '*':

					while (input[j] == " ")
					{
						j--;
					}
					a = stold(input[j]);
					whatOperation.insert(0, input[j]);
					whatOperation.insert(input[j].size(), " *");
					while (input[k] == " ")
					{
						k++;
					}
					b = stold(input[k]);
					whatOperation.insert(input[j].size() + 2, " " + input[k]);
					sequence[order] = whatOperation; //тупо добавляем в список операций строку, что делали
					order++;
					whatOperation.clear();
					localResult += round(a * b);

					input.erase(input.begin() + j, input.begin() + k + 1); //удаляем, что посчитали 

					break;

				case '^':

					while (input[j] == " ")
					{
						j--;
					}
					a = stold(input[j]);
					whatOperation.insert(0, input[j]);
					whatOperation.insert(input[j].size(), " ^");
					while (input[k] == " ")
					{
						k++;
					}
					b = stold(input[k]);
					whatOperation.insert(input[j].size() + 2, " " + input[k]);
					sequence[order] = whatOperation; //тупо добавляем в список операций строку, что делали
					order++;
					whatOperation.clear();
					localResult += round(pow(a, b));
					input.erase(input.begin() + j, input.begin() + k + 1); //удаляем, что посчитали 

					break;

				case '/':

					while (input[j] == " ")
					{
						j--;
					}
					a = stold(input[j]);
					whatOperation.insert(0, input[j]);
					whatOperation.insert(input[j].size(), " /");
					while (input[k] == " ")
					{
						k++;
					}
					b = stold(input[k]);
					whatOperation.insert(input[j].size() + 2, " " + input[k]);
					sequence[order] = whatOperation; //тупо добавляем в список операций строку, что делали
					order++;
					whatOperation.clear();
					localResult += round(a / b);
					input.erase(input.begin() + j, input.begin() + k + 1); //удаляем, что посчитали 

					break;

				case '%':

					long long _a, _b;
					while (input[j] == " ")
					{
						j--;
					}
					_a = stol(input[j]);
					whatOperation.insert(0, input[j]);
					whatOperation.insert(input[j].size(), " %");
					while (input[k] == " ")
					{
						k++;
					}
					_b = stol(input[k]);
					whatOperation.insert(input[j].size() + 2, " " + input[k]);
					sequence[order] = whatOperation; //тупо добавляем в список операций строку, что делали
					order++;
					whatOperation.clear();
					localResult += round(_a % _b);
					input.erase(input.begin() + j, input.begin() + k + 1); //удаляем, что посчитали 

					break;
				}
				input.emplace(input.begin() + j, to_string(localResult)); //вместо удалённой операции вставляем результат
			}

		}

		for (size_t i = 0; i < input.size(); i++)
		{ //надо рассмотреть отдельный случай, когда минус стоит в начале. вся фигня со скобками нестрашна, т.к. рекурсивно это всё равно сводится к случаю "минус в начале"
			long double localResult = 0;
			if ((input[i] == "-") or (input[i] == "+"))
			{
				long double a = 0.0;
				long double b = 0.0;
				size_t j = i - 1, k = i + 1;

				switch (input[i][0])
				{
				case '+':
					while (input[j] == " ")
					{
						j--;
					}
					a = stol(input[j]);
					whatOperation.insert(0, input[j]);
					whatOperation.insert(input[j].size(), " +");
					while (input[k] == " ")
					{
						k++;
					}
					b = stol(input[k]);
					whatOperation.insert(input[j].size() + 2, " " + input[k]);
					sequence[order] = whatOperation; //тупо добавляем в список операций строку, что делали
					order++;
					whatOperation.clear();
					localResult += round(a + b);
					input.erase(input.begin() + j, input.begin() + k + 1); //удаляем, что посчитали 

					break;

				case '-':
					if (i == 0)
					{
						//тот самый противный случай
						while (input[k] == " ")
						{
							k++;
						}
						input[k].insert(0, "-"); //явный костыль

						break;
					}
					else
					{
						while (input[j] == " ")
						{
							j--;
						}
						a = stol(input[j]);
						whatOperation.insert(0, input[j]);
						whatOperation.insert(input[j].size(), " -");
						while (input[k] == " ")
						{
							k++;
						}
						b = stol(input[k]);
						whatOperation.insert(input[j].size() + 2, " " + input[k]);
						sequence[order] = whatOperation; //тупо добавляем в список операций строку, что делали
						order++;
						whatOperation.clear();
						localResult += round(a - b);
						input.erase(input.begin() + j, input.begin() + k + 1); //удаляем, что посчитали 

						break;

					}
				}
				input.insert(input.begin() + j, to_string(localResult)); //вместо удалённой операции вставляем результат
			}
		}

		if (input.size() == 1)
		{
			result = round(stol(input[0]));
			return;
		}
		else if (input.size() == 2)
		{
			bool isEmpty = false;
			size_t notEmpty = 0;
			for (size_t i = 0; i < input.size(); i++)
			{
				if (input[i] == "")
				{
					isEmpty = true;
					notEmpty = input.size() - i - 1;
					result = round(stol(input[notEmpty]));
					return;
				}
			}
		}
		else if (input.size() == 0)
		{
			cerr << "Error in calculations! Result writing causes a problem!\n";
			exit(1);
		}

	}

}

bool const FileReadCalc::check(char symb)
{
	if ((static_cast<int>(symb) > 47) and (static_cast<int>(symb) < 58) or symb == '.' or symb == '(' or symb == ')')
	{
		return true;
	}
	return false;
}
