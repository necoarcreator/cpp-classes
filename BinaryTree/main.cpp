#include <iostream>
#include "CBinaryTree.h"

using namespace std;

int main()
{
	CBinaryTree<int> A{ 5 };
	A.addKnot(4);
	A.addKnot(6);
	A.addKnot(9);
	A.addKnot(7);
	A.addKnot(1);

	A.addKnot(2);
	A.addKnot(3);
	A.addKnot(15);
	A.addKnot(8);
	cout << A.delKnot(8) << endl;
	cout << "________________________________________" << endl;
	A.print(A.getroot());
	return 0;
}
