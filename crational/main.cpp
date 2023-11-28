#include "Crational.h"
#include <iostream>

using namespace std;

int main()
{
	float aa = 1.500;
	Crational a{ aa };
	float c = 2.5;
	float b = a * a;
	cout << a.getnum() << "/" << a.getdenom() << endl;
	cout << b;

	return 0;
} 