#include "Queue.h"
#include <vector>
#include <iostream>

using namespace std;

Queue::Queue(vector<int> init) : Datastruct(init) {}

Queue::Queue() : Datastruct() {}

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

void Queue::print()
{
    for (int i = 0; i < size; i++)
    {
        cout << A[i] << " ";
    }
    cout << endl;
}