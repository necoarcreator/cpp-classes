#include <iostream>
#include <vector>
#include <string>
#include <conio.h>
#include <windows.h>
#include <fstream>
#include <map>
#include <algorithm>
#include <locale>
#include <functional>

using namespace std;

int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	ifstream fs("tolstoj_lew_nikolaewich-text_0070.txt");


	map<const char, long long int> tree;
	
	string slovar{ "абвгдеёжзийклмнопрстуфхцчшщъыьэюя" };

	for (int i = 0; i < size(slovar); i++)
	{
		tree[slovar[i]] = 0;
	}
	
	while (fs)
	{
		
		string line;
		fs >> line;
		
		for (int i = 0; i < size(line); i++)
		{
			if ((static_cast<int>(line[i]) > -65) and (static_cast<int>(line[i]) < -32))
			{
				line[i] = line[i] + 'а' - 'А'; //аналог tolower для кириллицы
			}
		}
		for (int j = 0; j < size(line); j++) //смотрим всю строку
		{
			for (const auto letter : slovar) //проходимся по алфавиту
			{
				if (line[j] == letter) //если они равны
				{
					tree[letter] += 1;
				}
			}
		}

	}
	
	
	//cout << static_cast<int>('А') << endl;
	//cout << static_cast<int>('Я') << endl;

	for (auto it = tree.cbegin(); it != tree.cend(); ++it)
	{
		cout << "Буква " << it->first << " встретилась " << it->second << " раз.\n";
	}
	
	map <int, string> top;
	
	cout << "_____________________________________\n";
	cout << "Топ встречаемых символов:\n";

	for (int i = 1; i < 11; i++)
	{	
		long long int max = -1;
		char iter;
		for (auto it = tree.cbegin(); it != tree.cend(); ++it)
		{
			if (it->second > max)
			{
				iter = it->first;
				max = it->second;
			}
		}
		top[i] = iter;
		tree.erase(tree.find(iter));
	}
	
	for (auto it = top.cbegin(); it != top.cend(); ++it)
	{
		cout << "Топ - " << it->first << " буква по частоте - " << it->second << ".\n";
	}
	
	
	fs.close();

	return 0;
}