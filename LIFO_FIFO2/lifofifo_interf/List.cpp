#include "List.h"
#include <vector>
#include <iostream>
using namespace std;
template <typename T>

List<T>::List(vector<T> init) : A(init), size(A.size())
{
	for (int i = 1; i < size; i++)
	{
		Aptr.push_back(&A[i]);
	}

	Aptr.push_back(nullptr);
	print();
	printPtr();
}

template <typename T>

List<T>::List() : A{}, size(0), Aptr{}
{
	cout << "Warning! The list is empty.\n";
}

template <typename T>

void List<T>::push(int num)
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

template <typename T>

bool List<T>::isEmpty()
{
	if (size == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

template <typename T>

int List<T>::pop()
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

template <typename T>

int List<T>::search(T x)
{
	int i = 0;
	while (Aptr[i] != nullptr && A[i] != x)
	{
		i++;
	}
	return i;
}

template <typename T>

int List<T>::popFront()
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

template <typename T>

List<T> List<T>::rev()
{
	if (not isEmpty())
	{
		for (int i = 0; i < int(size / 2); i++)
		{
			swap(A[i], A[size - i - 1]);
		}

		Aptr.clear();
		Aptr.shrink_to_fit();
		for (int i = 1; i < size; i++)
		{
			Aptr.push_back(&A[i]);
		}

		Aptr.push_back(nullptr);

		return A;
	}
	else
	{
		cout << "Nothing to swap!\n";
		
	}
}

template <typename T>

void List<T>::print()
{
	for (int i = 0; i < size; i++)
	{
		cout << A[i] << " ";
	}
	cout << endl;
}

template <typename T>

void List<T>::printPtr()
{
	for (int i = 1; i < size; i++)
	{
		cout << *Aptr[i - 1] << " __ ";
	}
	cout << endl;
}

template <typename T>

vector<T> List<T>::getA()
{
	return A;
}

template <typename T>

int List<T>::getsize()
{
	return size;
}