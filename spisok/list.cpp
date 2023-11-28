#include "list.h"
#include <vector>
#include <iostream>
using namespace std;


List::List(vector<int> init) : A(init), size(A.size()) 
{
	for (int i = 1; i < size; i++)
	{
		Aptr.push_back(&A[i]);
	}

	Aptr.push_back(nullptr);
	print();
	printPtr();
}

List::List() : A{}, size(0), Aptr{}
{
	cout << "Warning! The list is empty.\n";
}

void List::push(int num)
{
	A.push_back(num);
	Aptr.clear();
	Aptr.shrink_to_fit();

	size += 1;
	for (int i = 1; i < size; i++)
	{
		Aptr.push_back(&A[i]);
	}

	Aptr.push_back(nullptr);
	print();
	printPtr();

	return;
}

int List::pop()
{
	if (size == 0)
	{
		cout << "Error: no data to delete!\n";
		return 0;
	}

	int num = A[size - 1];
	Aptr.erase(Aptr.end() - 1);
	A.erase(A.end() - 1);
	Aptr[size - 2] = nullptr;
	size -= 1;
	printPtr();
	print();
	return num;
}

int List::pop_front()
{
	if (size == 0)
	{
		cout << "Error: no data to delete!\n";
		return 0;
	}
	int num = A[0];
	A.erase(A.begin());
	Aptr.erase(Aptr.end() - 1);
	Aptr[size - 2] = nullptr;
	size -= 1;
	printPtr();
	print();
	return num;
}
void List::print()
{
	for (int i = 0; i < size; i++)
	{
		cout << A[i] << " ";
	}
	cout << endl;
}

void List::printPtr()
{
	for (int i = 1; i < size; i++)
	{
		cout << *Aptr[i - 1] << " __ ";
	}
	cout << endl;
}