#include <iostream>
#include "CBinaryTree.h"

using namespace std;

int main()
{
	CBinaryTree<double> A{ 5.7 };
	A.addKnot(4.5);
	A.addKnot(6.1);
	A.addKnot(9.4);
	A.addKnot(7.9);
	A.addKnot(1.2);
	
	A.addKnot(2.4);
	A.addKnot(3.6);
	A.addKnot(15.1);
	cout << A.delKnot(3.6) << endl;
	A.addKnot(8.9);
	cout << "________________________________________" << endl;
	A.print(A.getroot());
	return 0;
}
