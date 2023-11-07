#ifndef _QE_H
#define _QE_H
#include <vector>
using namespace std;
class Queue
{
    vector<int> A;
    int size;

public:
    Queue();

    Queue(vector<int> inpt);

    void push(int arg);

    int pop();

    void print();
};

#endif