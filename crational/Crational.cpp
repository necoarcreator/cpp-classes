#include "Crational.h"
#include <math.h>
#include <iostream>
#include <limits>
#include <string>
#include "CrationalExeption.h"
using namespace std;


Crational::Crational(float ins) : denom{1}, num{int(ins)}
{
	
	if (ins > 0)
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
	}
	else
	{
		float modins = -ins;
		modins -= num;
		float dr = num;
		num = 0;
		int count = 0;

		while ((modins - int(modins) > 1e-4) and (count < 8))
		{
			modins *= 10;
			num = num * 10 + int(modins);
			modins -= int(modins);

			denom *= 10;
			count += 1;
		}
		num = num + pow(10, count) * dr;
		num *= -1;
	}
};

Crational::Crational(int chis, int znam) : num{ chis }, denom{ znam }
{
	if (denom == 0)
		{
			throw CrationalExeption("Error! Division by zero!");
		}
		else if (abs(chis / znam) >= FLT_MAX)
		{
			throw - 1.0;
		}
}

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

float Crational::operator- (Crational b)
{	

	float k = (num * b.getdenom() - b.getnum() * denom);
	k /= denom * b.getdenom();
	return k;
}

float Crational::operator- (float b)
{
	Crational bcr{ b };
	float k = (num * bcr.getdenom() - bcr.getnum() * denom);
	k /= denom * bcr.getdenom();
	return k;
}

float Crational::operator- (int b)
{
	float k = (num - denom * b);
	k /= denom;
	return k;
}

float Crational::operator/ (int a)
{
	return denom * a;
}

float Crational::operator/ (float a)
{
	Crational acr{ a };
	float k = (denom * acr.getnum());
	k /= num * acr.getdenom();
	return k;
}

float Crational::operator/ (Crational a)
{
	float k = (denom * a.getnum());
	k /= (num * a.getdenom());
	return k;
}

Crational::~Crational()
{
	num = 0;
	denom = 0;

}