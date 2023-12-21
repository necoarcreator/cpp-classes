#include "Stack.h"
#include <iostream>
#include <vector>

using namespace std;

Stack::Stack() : size(0), A{}
{
    cout << "Warning: the stack is empty!\n";
}

Stack::Stack(vector<int> inpt) : A(inpt), size(A.size())
{
    print();
}

void Stack::push(int arg)
{
    A.push_back(arg);
    size += 1;
    print();
    return;
}

int Stack::pop()
{
    try
    {
        if (size == 0)
        {
            throw 1;
        }
    }
    catch (const int exp)
    {
        cout << "Error: " << exp << ". No data to delete!\n";
        return 0;
    }
    int arg = A[size - 1];
    A.pop_back();
    size -= 1;
    print();
    return arg;
}

bool Stack::isEmpty()
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

Stack Stack::rev()
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

void Stack::print()
{
    for (int i = 0; i < size; i++)
    {
        cout << A[i] << " ";
    }
    cout << endl;
}

vector<int> Stack::getA()
{
    return A;
}

int Stack::getsize()
{
    return size;
}