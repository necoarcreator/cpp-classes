#include "WeightGraph.h"
#include <iostream>
#include <list>
using namespace std;
WeightGraph::WeightGraph()
{
    pknots = new list<Knot>;
}
WeightGraph::WeightGraph(int numberKnots)
{
    pknots = new list<Knot>;
    for (int i = 0; i < numberKnots; i++)
    {
        this->addKnot(i);
    }
}
WeightGraph::WeightGraph(WeightGraph const& b) 
{
    pknots = new list<Knot>(*(b.pknots));
};
void WeightGraph::addKnot(int index)
{
	pknots->emplace_back(index);
}
void WeightGraph::addKnot(Knot index)
{
    pknots->emplace_back(index);
}

void WeightGraph::addEdge(int source, int dest, int weight)
{
    for (Knot& u : *pknots) 
    {
            if (*(u.porder) == source)
            {
                
                u.pedges->emplace_back(dest, weight);
            }
            else if (*(u.porder) == dest)
            {
                u.pedges->emplace_back(source, weight);
            }
        
    }
}
void WeightGraph::deleteEdge(int _Knot, int _Edge)
{
    for (Knot& u : *pknots)
    {
        if (*(u.porder) == _Knot)
        {
            auto iter = u.pedges->begin();
            advance(iter, _Edge);
            u.pedges->erase(iter);
        }
    }
}

void WeightGraph::printGraph()
{
    for (Knot& u : *pknots)
    {
        cout << "Knot number " << *u.porder << " has connection with: " << endl;

        for (Edge& v : *(u.pedges))
        {
            cout << "knot number " << *v.pdest << " with weight " << *v.pweight << endl;
        }

        cout << endl;
    }
}

Knot WeightGraph::getKnots(int index)
{
    auto iter = pknots->begin();
    advance(iter, index);
    return (*iter);
}

Knot WeightGraph::getKnots(int* porder)
{
    auto iter = pknots->begin();
    for (unsigned int i = 0; i < pknots->size(); i++)
    {
        if (*(getKnots(i).porder) == *porder)
        {
            advance(iter, i);
            return (*iter);
        }
    } 
    cerr << "Error! Knot was not found";
    exit(1);
}

Edge WeightGraph::getEdges(int knotIndex, int edgeIndex)
{
    auto iter1 = pknots->begin();
    advance(iter1, knotIndex);
    auto iter2 = iter1->pedges->begin();
    advance(iter2, edgeIndex);
    return (*iter2);
}


int WeightGraph::knotsSize()
{
    return pknots->size();

}

list<Knot>* WeightGraph::copyKnots()
{
    return pknots;
}


WeightGraph WeightGraph::operator+(WeightGraph b)
{
    int n(knotsSize());
    int m(b.knotsSize());


    if (n >= m)
    {
        WeightGraph* c = new WeightGraph;
        for (int k = 0; k < n; k++)
        {
            (*c).addKnot(getKnots(k));
        }

        for (int k = 0; k < n; k++)
        {
            for (unsigned int q = 0; q < getKnots(k).pedges->size(); q++)
            {
                bool unique = true;
                for (unsigned int p = 0; p < (*c).getKnots(k).pedges->size(); p++)
                {
                    int n1 = *(getEdges(k, q)).pdest;
                    int n2 = *((*c).getEdges(k, p)).pdest;

                    if ((n1) == (n2))
                    {
                        unique = false;
                        break;
                    }
                }
                if (unique)
                {
                    int d = *(getKnots(k).porder);
                    int nn = *(getEdges(k, q)).pdest;
                    int ww = *(getEdges(k, q)).pweight;
                    (*c).addEdge(d, nn, ww);
                }

            }
        }

        for (int k = 0; k < m; k++)
        {
            int mm = *(b.getKnots(k).porder);
            bool unique = true;
            for (int q = 0; q < n; q++)
            {
                int nn = *(getKnots(q).porder);
                if (nn == mm)
                {
                    unique = false;
                    break;
                }
            }
            if (unique)
            {
                (*c).addKnot(*(b.getKnots(k).porder));
            }

        }

        int s((*c).knotsSize());
        
            for (int k = 0; k < m; k++)
            {
                for (int j = 0; j < s; j++)
                {
                    if (*(b.getKnots(k).porder) == *((*c).getKnots(j).porder))
                    {
                        for (unsigned int q = 0; q < b.getKnots(k).pedges->size(); q++)
                        {
                            bool unique = true;

                            for (unsigned int p = 0; p < (*c).getKnots(j).pedges->size(); p++)
                            {
                                int n1 = *(b.getEdges(k, q)).pdest;
                                int n2 = *((*c).getEdges(j, p)).pdest;

                                if ((n1) == (n2))
                                {
                                    unique = false;
                                    break;
                                }
                            }
                            if (unique)
                            {
                                int d = *((*c).getKnots(j).porder);
                                int nn = *(b.getEdges(k, q)).pdest;
                                int ww = *(b.getEdges(k, q)).pweight;
                                (*c).addEdge(d, nn, ww);
                            }
                        }
                    }


                }
            }
            WeightGraph d = *c;
            delete c;
        return d;
    }
    else
    {
        return b + (*this);
        
    }
}

WeightGraph WeightGraph::operator-(WeightGraph b)
{
    int n(knotsSize());
    int m(b.knotsSize());

    WeightGraph c;
    if (n >= m)
    {

        for (int k = 0; k < n; k++)
        {
            c.addKnot(k);
            for (unsigned int q = 0; q < getKnots(k).pedges->size(); q++)
            {
                auto iter = getKnots(k).pedges->begin();
                advance(iter, q);
                int nn = *(*iter).pdest;
                int ww = *(*iter).pweight;
                c.addEdge(k, nn, ww);
            }
        }

        for (int k = 0; k < n; k++)
        {
            for (unsigned int q = 0; q < c.getKnots(k).pedges->size(); q++)
            {
                bool unique = true;
                for (unsigned int p = 0; p < c.getKnots(k).pedges->size(); p++)
                {
                    auto iterc1 = c.getKnots(k).pedges->begin();
                    auto iterc2 = c.getKnots(k).pedges->begin();
                    advance(iterc1, q);
                    advance(iterc2, p);
                    if (((*iterc2).pdest) == ((*iterc2).pdest) && (q != p))
                    {
                        unique = false;
                        break;
                    }
                }
                if (not unique)
                {
                    c.deleteEdge(k, q);
                }
            }
        }
    }
    else
    {

        return b - (*this);

    }
    cerr << "Unexpected error!";
    exit(1);
}
WeightGraph& WeightGraph::operator=(WeightGraph b)
{
    swap(*pknots, *(b.copyKnots()));

    return *this;
}



WeightGraph::~WeightGraph()
{
    delete pknots;
    pknots = nullptr;
}