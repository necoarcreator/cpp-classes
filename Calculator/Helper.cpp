#include "Helper.h"

#include <iostream> 
#include <fstream>
#include <string>
#include <vector>

using namespace std;
Helper::Helper()
{

	ifstream file("commandList.txt");
	if (file.is_open())
	{
		string line;

		while (getline(file, line))
		{
			output.push_back(line);
			cout << line << endl;
		}
	}
	else
	{
		cerr << "Can't open file \"commandList.txt\" for reading" << endl;
		exit(1);
	}
	file.close();

}

Helper::~Helper()
{
	output.clear();

}

void const Helper::addCommand(string input)
{
	ofstream out("commandList.txt");
	if (out.is_open())
	{
		out << input;
		output.push_back(input);
		out << "\n";
	}
	else
	{
		cerr << "Can't open file \"commandList.txt\" for writing\n";
	}
	out.close();
}
