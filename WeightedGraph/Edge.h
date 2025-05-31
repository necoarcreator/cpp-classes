#ifndef _ED_
#define _ED_

#include <vector>
using namespace std;
struct Edge
{
	int* pdest;
	int* pweight;

	Edge();
	Edge(int _dest, int _weight);
	Edge(Edge const& b);
	Edge& operator=(Edge b);
	~Edge();
};

#endif