#include "Queue.h"
#include <iostream>
#include <vector>

using namespace std;



Queue::Queue() : size(0)
{
    A.push_back(0);
    cout << "Warning: the stack is empty!\n";
}

Queue::Queue(vector<int> inpt) : A(inpt), size(inpt.size())
{}

void Queue::push(int arg)
{
    A.push_back(arg);
    size += 1;
    print();
    return;
}

int Queue::pop()
{
    if (size == 0)
    {
        cout << "Error: no data to delete!\n";
        return 0;
    }
    int arg = A[0];
    A.erase(A.begin());
    size -= 1;
    print();
    return arg;
}

void Queue::print()
{
    for (int i = 0; i < size; i++)
    {
        cout << A[i] << " ";
    }
    cout << endl;
}