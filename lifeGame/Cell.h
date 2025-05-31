#ifndef _CELL_
#define _CELL_

struct Cell
{
	bool state, nextState;
	
	Cell(bool _state);
	Cell();
	~Cell();


};


#endif