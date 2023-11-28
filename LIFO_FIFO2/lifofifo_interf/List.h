#ifndef _LST_

#define _LST_

#include "interface.h"
#include <vector>

using namespace std;
template <typename T>

class List : public Interface
{
private:
	vector<T> A;
	int size;
	vector<T*> Aptr;
public:

	List(vector<T> init);

	List();

	void push(int num) override final;

	int pop() override final;

	bool isEmpty() override final;

	List rev();

	int search(T x);

	int popFront();

	vector<T> getA();

	int getsize();

	void print();

	void printPtr();

};
#endif 