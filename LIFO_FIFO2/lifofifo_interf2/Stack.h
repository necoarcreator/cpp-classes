#ifndef _STK_H
#define _STK_H

#include "Datastruct.h"
#include <vector>

using namespace std;

class Stack : public Datastruct
{
public:
	Stack(vector<int> init);

	Stack();

	void push(int num) override;

	int pop() override;

	void print();
};

#endif