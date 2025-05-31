#include "Cell.h"

Cell::Cell(bool _state) : state(_state), nextState(_state) {};
Cell::Cell() : state(false), nextState(false) {};
Cell::~Cell()
{
	state = false;
	nextState = false;
}