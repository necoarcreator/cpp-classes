#include "list.h"
#include <iostream>
#include <vector>

using namespace std;

int main()
{
    vector<int> a{ 1,2,3,4,5,6,7 };
    List L(a);

    L.pop_front();

}

