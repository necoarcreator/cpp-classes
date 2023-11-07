#include "Stack.h"
#include <vector>
#include <iostream>

using namespace std;

Stack::Stack(vector<int> init) : Datastruct(init) {}

Stack::Stack() : Datastruct() {}

void Stack::push(int num)
{
    A.push_back(num);
    size += 1;
    print();
    return;
}

int Stack::pop()
{
    if (size == 0)
    {
        cout << "Error: no data to delete!\n";
        return 0;
    }
    int arg = A[size - 1];
    A.pop_back();
    size -= 1;
    print();
    return arg;
}

void Stack::print()
{
    for (int i = 0; i < size; i++)
    {
        cout << A[i] << " ";
    }
    cout << endl;
}