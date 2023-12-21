#ifndef _EXP

#define _EXP

#include <string>
using namespace std;

class CrationalExeption
{
private:
	const char* whatToCout;
public:
	CrationalExeption(const char divZero[]);

	CrationalExeption(float overFlow);
	
	const char* getError();
};

#endif 