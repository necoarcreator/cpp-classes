#include "List.h"
#include <iostream>
using namespace std;

#include "List.h"

List::List() : head(nullptr), tail(nullptr) {};

void List::push(int num) {
	Node* p = new Node(num);
	if (head == NULL) {
		head = p;
	}
	if (tail != NULL) {
		tail->next = p;
	}
	tail = p;
}

int List::pop() {
	int temp;
	try
	{
		if (tail == NULL) {
			throw 1;
		}
	}
	catch (const int exp)
	{
		cout << "Error: " << exp << ". No data to delete!\n";
		return 0;
	}
	if (head == tail) {
		temp = tail->val;
		delete tail;
		head = tail = NULL;
		return temp;
	}

	Node* p = head;
	while (p->next != tail) {
		p = p->next;
	}

	p->next = NULL;
	temp = tail->val;
	delete tail;
	tail = p;
	return temp;
}

void List::print() {
	if (head == NULL) {
		std::cout << "The list is empty!";
	}
	Node* p = head;
	while (p) {
		std::cout << p->val << " ";
		p = p->next;
	}
	std::cout << "\n";
}

bool List::isEmpty()
{
	if (head == NULL)
	{
		return true;
	}
	return false;
}
Node* List::getTail()
{
	return tail;
}
/*
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


bool List::isEmpty()
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


int List::search(int x)
{
	int i = 0;
	while (Aptr[i] != nullptr && A[i] != x)
	{
		i++;
	}
	return i;
}


int List::popFront()
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


List List::rev()
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


vector<int> List::getA()
{
	return A;
}


int List::getsize()
{
	return size;
}

*/