#ifndef _WG_
#define _WG_

#include <list>
#include "Edge.h"
#include "Knot.h"

using namespace std;

class WeightGraph
{
	list<Knot>* pknots;
	
public:

	WeightGraph();

	WeightGraph(int numberKnots);

	WeightGraph(WeightGraph const& b);

    void addKnot(int index);

	void addKnot(Knot index);

	void addEdge(int source, int dest, int weight);

	void deleteEdge(int _Knot, int _Edge);

	void printGraph();

	Knot getKnots(int index);

	Knot getKnots(int* porder);

	Edge getEdges(int knotIndex, int edgeIndex);

	list<Knot>* copyKnots();

	int knotsSize();

	WeightGraph operator+(WeightGraph b);

	WeightGraph operator-(WeightGraph b);

	WeightGraph& operator=(WeightGraph b);

	~WeightGraph();

};

#endif