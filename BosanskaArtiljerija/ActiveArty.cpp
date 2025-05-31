#include "ActiveArty.h"

ActiveArty::ActiveArty(double _range, double _temperature, string _type, double _caliber, size_t _numInBattery, size_t _numShells,
	string _typeshells, string _typefuze, string _name) : range(_range), temperature(_temperature), typesCalibers(make_pair(_type, _caliber)),
	numInBattery(_numInBattery), numShells(_numShells), typeshells(_typeshells), typefuze(_typefuze), name(_name) {};

ActiveArty::~ActiveArty()
{
	return;
}

ActiveArty::ActiveArty() : range(0), temperature(0), typesCalibers(make_pair("", 0.)),
numInBattery(0), numShells(0), typeshells(""), typefuze(""), name ("") {};