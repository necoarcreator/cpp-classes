#include "Queue.h"
#include <vector>
#include <iostream>

using namespace std;

Queue::Queue(vector<int> init) : A(init), size(A.size())
{
    print();
}

Queue::Queue() : size(0), A{}
{
    cout << "Warning: the queue is empty!\n";
}

void Queue::push(int num)
{
    A.push_back(num);
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
    int num = A[0];
    A.erase(A.begin());
    size -= 1;
    print();
    return num;
}

bool Queue::isEmpty()
{
    if (size == 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

Queue Queue::rev()
{
    if (not isEmpty())
    {
        for (int i = 0; i < int(size / 2); i++)
        {
            swap(A[i], A[size - i - 1]);
        }

        return A;
    }
    else
    {
        cout << "Nothing to swap!\n";
    }
}

void Queue::print()
{
    for (int i = 0; i < size; i++)
    {
        cout << A[i] << " ";
    }
    cout << endl;
}

vector<int> Queue::getA()
{
    return A;
}

int Queue::getsize()
{
    return size;
}