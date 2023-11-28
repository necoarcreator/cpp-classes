#include "Crational.h"
#include <math.h>
using namespace std;


Crational::Crational(float ins) : denom{1}, num{int(ins)}
{
	ins -= num;
	float dr = num;
	num = 0;
	int count = 0;

	while ((ins - int(ins) > 1e-4) and (count < 8))
	{	
		ins *= 10;
		num = num * 10 + int(ins);
		ins -= int(ins);
		
		denom *= 10;
		count += 1;
	}
	num = num + pow(10, count) * dr;
};

Crational::Crational(int chis, int znam) : num{ chis }, denom{ znam } {};

int Crational::getnum()
{
	return num;
}

int Crational::getdenom()
{
	return denom;
}

float Crational::operator+ (Crational b)
{
	float k = (num * b.getdenom() + b.getnum() * denom);
	k /= denom * b.getdenom();
	return k;
}

float Crational::operator+ (float b)
{
	Crational bcr{ b };
	float k = (num * bcr.getdenom() + bcr.getnum() * denom);
	k /= denom * bcr.getdenom();
	return k;
}

float Crational::operator+ (int b)
{
	float k = (num + denom * b);
	k /= denom;
	return k;
}

float Crational::operator* (int a)
{
	return num * a;
}

float Crational::operator* (float a)
{
	Crational acr{ a };
	float k = (num * acr.getnum());
	k /= denom * acr.getdenom();
	return k;
}

float Crational::operator* (Crational a)
{
	float k = (num * a.getnum());
	k /= (denom * a.getdenom());
	return k;
}