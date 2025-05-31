#include "Knot.h"
#include <iostream>
Knot::Knot()
{
	porder = new int;
	pedges = new list<Edge>;
};

Knot::Knot(int i)
{
	porder = new int;
	*porder = i;
	pedges = new list<Edge>;

};
Knot::Knot(Knot const & b) 
{
	porder = new int(*(b.porder));
	pedges = new list<Edge>(*(b.pedges));
};
Knot& Knot::operator=(Knot b)
{
	swap(*pedges, *(b.pedges));
	swap(*porder, *(b.porder));
	return *this;
}
Knot::~Knot()
{
	
	delete porder;
	porder = nullptr;
	
	delete pedges;
	pedges = nullptr;
	

};