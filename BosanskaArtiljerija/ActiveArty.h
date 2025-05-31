#ifndef _AA_
#define _AA_
#include <string>

using namespace std;

struct ActiveArty
{
	double range;
	double temperature;
	pair <string, double> typesCalibers;
	size_t numInBattery;
	size_t numShells;
	string typeshells;
	string typefuze;
	string name;

	ActiveArty(double _range, double _temperature, string _type, double _caliber, size_t _numInBattery, size_t _numShells,
		string _typeshells, string _typefuze, string _name);
	ActiveArty();
	~ActiveArty();
};

#endif