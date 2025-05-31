#include "Edge.h"

Edge::Edge() 
{
	pdest = new int;
	pweight = new int;
};
Edge::Edge(int _dest, int _weight) 
{
	pdest = new int;
	*pdest = _dest;
	pweight = new int;
	*pweight = _weight;

}
Edge::Edge(Edge const& b) 
{
	pdest = new int(*(b.pdest));
	pweight = new int(*(b.pweight));
};
Edge& Edge::operator=(Edge b)
{
	swap(*pweight, *(b.pweight));
	swap(*pdest, *(b.pdest));
	return *this;
}
Edge::~Edge()
{
	
	delete pdest;
	pdest = nullptr;
	delete pweight;
	pweight = nullptr;
}