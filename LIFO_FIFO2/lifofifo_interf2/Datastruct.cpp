
#include "Datastruct.h"
#include <vector>

using namespace std;

Datastruct::Datastruct(vector<int> init) : A(init), size(init.size()) {}
Datastruct::Datastruct() : A({}), size(0) {}