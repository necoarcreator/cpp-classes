#include "StandartCalc.h"

#include <iostream>
#include<array>
#include <string>
#include <vector>
#include <map>
#include <math.h>
using namespace std;

StandartCalc::StandartCalc() : listOfCommands("+-*/%^"), argc(0), order(1), numRequests(0)
{
}
StandartCalc::StandartCalc(int _argc, vector<string> input, bool isWritingLog) : listOfCommands("+-*/%^"), argc(_argc), order(1)
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
	{//
		originalString.push_back({});
		originalString[i - 2 - movingOfIterator].push_back("");
		size_t p = 0; //счетчик длины числа
		size_t j = 0; //счетчик числа элементов в строке ввода
		for (size_t k = 0; k < input[i].size(); k++)
		{//

			if (check(input[i][k]))
			{
				p++;

			}

			if ((input[i][k] == ' ') and (check(input[i][k - 1])))
			{
				string sub = input[i].substr(k - p, p);

				originalString[i - 2 - movingOfIterator][j] += sub;
				originalString[i - 2 - movingOfIterator].push_back("");
				j++;
				p = 0;
			}
			else if (not check(input[i][k]))
			{

				for (auto x : listOfCommands)
				{
					if (x == input[i][k])
					{
						if (check(input[i][k - 1]))
						{
							string sub = input[i].substr(k - p, p);

							originalString[i - 2 - movingOfIterator][j] += sub;
							originalString[i - 2 - movingOfIterator].push_back("");
							j++;
							p = 0;
						}

						originalString[i - 2 - movingOfIterator][j] += input[i][k]; 
						originalString[i - 2 - movingOfIterator].push_back("");
						j++;

						if (check(input[i][k + 1]))
						{
							p++;
						}
						p = 0;
					}
				}

			}

			if (k == input[i].size() - 1)
			{
				string sub = input[i].substr(k - p + 1, p);

				originalString[i - 2 - movingOfIterator][j] += sub;
			}

		}
	}

	numRequests = input.size() - 2 - movingOfIterator;

}
void const StandartCalc::printSequence()
{
	cout << "********************************" << endl;
	for (auto it = sequence.begin(); it != sequence.end(); ++it)
	{
		cout << "operation number " << it->first << ": " << it->second << endl;
		cout << "********************************" << endl;
	}
}
StandartCalc::~StandartCalc()
{
	originalString.clear();
	sequence.clear();
}

bool const StandartCalc::check(char symb)
{
	if ((static_cast<int>(symb) > 47) and (static_cast<int>(symb) < 58) or symb == '.' or symb == '(' or symb == ')')
	{
		return true;
	}
	return false;
}

void const StandartCalc::startCount()
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
void const StandartCalc::proceed(vector<string> input)
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




