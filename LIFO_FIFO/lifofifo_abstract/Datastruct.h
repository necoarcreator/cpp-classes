#ifndef _DTST_H
#define _DTST_H
#include "Struct.h"
#include <vector>

using namespace std;

class Datastruct : public Struct
{
protected:

	vector <int> A;
	int size;

public:

	Datastruct(vector<int> init);
	Datastruct();
};

#endif