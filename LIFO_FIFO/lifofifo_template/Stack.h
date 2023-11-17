#ifndef _ST_H
#define _ST_H
                                                                                                                                                        //#define A B
#include <iostream>
#include <vector>

using namespace std;

class Stack
{
    vector<int> A;
    int size;

public:
    Stack();

    Stack(int num, ...);

    Stack(vector<int> inpt);

    void push(int arg);

    int pop();

    void print();

};

#endif
