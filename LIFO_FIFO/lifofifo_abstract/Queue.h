#ifndef _QE_H
#define _QE_H

#include "Datastruct.h"
#include <vector>

using namespace std;

class Queue : public Datastruct
{
public:
	Queue(vector<int> init);

	Queue();

	void push(int num) override;

	int pop() override;

	void print();
};

#endif