#ifndef _LST_

#define _LST_
#include "Interface.h"
#include "Node.h"

class List : public Interface {
public:
	Node* head;
	Node* tail;

	List();
	virtual void push(int num);
	virtual int pop();
	virtual void print();
	bool isEmpty();
	Node* getTail();
};

#endif

/*
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

	void push(int num) override final;

	int pop() override final;

	bool isEmpty() override final;

	List rev();

	int search(int x);

	int popFront();

	vector<int> getA();

	int getsize();

	void print();

	void printPtr();

};
#endif 
*/