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
	
	string slovar{ "��������������������������������" };

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
				line[i] = line[i] + '�' - '�'; //������ tolower ��� ���������
			}
		}
		for (int j = 0; j < size(line); j++) //������� ��� ������
		{
			for (const auto letter : slovar) //���������� �� ��������
			{
				if (line[j] == letter) //���� ��� �����
				{
					tree[letter] += 1;
				}
			}
		}

	}
	
	
	//cout << static_cast<int>('�') << endl;
	//cout << static_cast<int>('�') << endl;

	for (auto it = tree.cbegin(); it != tree.cend(); ++it)
	{
		cout << "����� " << it->first << " ����������� " << it->second << " ���.\n";
	}
	
	map <int, string> top;
	
	cout << "_____________________________________\n";
	cout << "��� ����������� ��������:\n";

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
		cout << "��� - " << it->first << " ����� �� ������� - " << it->second << ".\n";
	}
	
	
	fs.close();

	return 0;
}