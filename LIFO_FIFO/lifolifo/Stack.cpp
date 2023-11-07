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
