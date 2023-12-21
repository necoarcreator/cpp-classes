#ifndef _TR_

#define _TR_

#include "Node.h"

using namespace std;
template <class T>

class CBinaryTree {



public:
	Node<T>* root;
	Node<T>* srchKnot;
	Node<T>** childKnot;
	Node<T>** fatherKnot;

	CBinaryTree();
	CBinaryTree(T num);
	void addKnot(T num);
	T delKnot(T num);
	T SrchPlc(T num);
	Node<T>* getroot();
	void print(Node<T>* root);
};









#endif

