#ifndef _ST_H
#define _ST_H

#include <iostream>
#include <vector>

using namespace std;

class Stack
{
protected:
    vector<int> A;
    int size;

public:
    Stack();

    Stack(vector<int> inpt);

    void push(int arg);

    int pop();

    void print();

};

#endif