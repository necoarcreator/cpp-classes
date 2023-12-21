#include "Crational.h"
#include "CrationalExeption.h"
#include <iostream>

using namespace std;

int main()
{
	float aa = 100000000000, bb = 0.000000000000000000000001;
	try {
		//Crational a{ 1, 0 };
		Crational a{ 1, 10 };
		Crational b{ -bb };
		float c = 2.5;
		Crational d{ aa / bb };
		cout << d.getnum() << "/" << d.getdenom() << endl;
	}
	catch (CrationalExeption &expn)
	{
		cerr << "An exeption occured: (" << expn.getError() << ")" << endl;
	}
	

	return 0;
} 