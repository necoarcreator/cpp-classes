#ifndef _QE_H
#define _QE_H

#include <iostream>
#include <vector>
#include "Stack.h"

using namespace std;

class Queue : public Stack
{

public:
    Queue();

    Queue(vector<int> inpt);

    int pop();

};

#endif
