#ifndef _KN_
#define _KN_

#include <list>
#include "Edge.h"
using namespace std;

struct Knot
{
	int* porder;

	list<Edge>* pedges;

	Knot();
	Knot(int i);
	Knot(Knot const& b);
	Knot& operator=(Knot b);
	~Knot();
};

#endif