#include <iostream>
#include "WeightGraph.h"

using namespace std;
int main()
{
	WeightGraph A, B;

    
    
    A.addKnot(1);
    
    A.addKnot(2);
    A.addKnot(3);
    A.addKnot(4);

    A.addEdge(1, 2, 5);
    A.addEdge(2, 3, 10);
    A.addEdge(1, 3, 7);
    A.addEdge(1, 4, 8);

    B.addKnot(2);
    B.addKnot(4);
    B.addKnot(5);
    B.addKnot(6);
    B.addKnot(7);

    B.addEdge(3, 4, 18);
    B.addEdge(4, 5, 20);
    B.addEdge(3, 6, 14);
    B.addEdge(1, 6, 16);
    B.addEdge(2, 7, 22);

    WeightGraph C = (A + B);
    

    C.printGraph();

}
