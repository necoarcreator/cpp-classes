#ifndef _LST_

#define _LST_

#include "interface.h"
#include <vector>

using namespace std;

class List : public Interface
{
private:
	vector<int> A;
	int size;
	vector<int*> Aptr;
public:

	List(vector<int> init);

	List();

	void push(int num) override;

	int pop() override;

	int pop_front();

	void print();

	void printPtr();
};
#endif 