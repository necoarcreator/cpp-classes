#include <iostream>
#include "crest.h"
#include <vector>

using namespace std; 

int main () {
	crest hop("input.txt");
	
	hop.solverCrestNotDiv("results1crest", ".csv", ',',
				10, 1, 2, 150);
	
	return 0;
}
