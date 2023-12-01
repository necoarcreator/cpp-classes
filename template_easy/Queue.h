#ifndef _QE_H
#define _QE_H

#include "interface.h"
#include <vector>

using namespace std;

class Queue : public Interface
{
private:
	vector<int> A;
	int size;

public:
	Queue(vector<int> init);

	Queue();

	void push(int num) override final;

	int pop() override final;

	bool isEmpty() override final;

	vector<int> getA();

	int getsize();

	Queue rev();

	void print();
};

#endif
