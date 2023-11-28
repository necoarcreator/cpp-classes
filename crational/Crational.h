#ifndef _RAT_
#define _RAT_
#include <minmax.h>
using namespace std;

//template <class T>

class Crational
{
private:

	int num, denom;

public:

	Crational(float ins);

	Crational(int chis, int znam);

	int getnum();

	int getdenom();

	float operator+(float b);

	float operator+(Crational b);

	float operator+(int b);

	float operator*(int a);

	float operator*(float a);

	float operator*(Crational a);

};




#endif