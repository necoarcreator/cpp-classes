#ifndef _STK_H
#define _STK_H

#include "interface.h"
#include <vector>

using namespace std;

class Stack : public Interface
{
private:
	vector<int> A;
	int size;
public:
	Stack(vector<int> init);

	Stack();

	void push(int num) override final;

	int pop() override final;

	bool isEmpty() override final;

	vector<int> getA();

	int getsize();

	Stack rev();

	void print();
};

#endif