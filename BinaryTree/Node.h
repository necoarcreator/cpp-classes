#ifndef _ND_
#define _ND_

using namespace std;

template <class T>
class Node
{
public:
	T val;
	/*T* left;
	T* right;
	left = new Node<T>();
	right = new Node<T>();*/
	Node<T>* left;
	Node<T>* right;

	Node();
	Node(T num);
};

#endif