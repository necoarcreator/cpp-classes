#include "Stack.h"
#include <iostream>
#include <vector>
#include <stdarg.h>
#include <limits.h>

using namespace std;

Stack::Stack() : size(0), A{}
{
    cout << "Warning: the stack is empty!\n";
}

Stack::Stack(int first = -INT_MIN, ...) : size(0) {
    if (first != -INT_MIN) {
        A.push_back(first);
        size += 1;
    }
    va_list factor;
    va_start(factor, first);
    while (*factor > 0) {
        A.push_back(va_arg(factor, int));
        size += 1;
    };
    va_end(factor);
}
/*
{
    va_list v1;
    int i;
    va_start (v1, num);

    for (i = 0; i < num; i++)
    {
        A.push_back(va_arg(v1, int));
        cout << "\n the number is" << A.back() << "\n";
    }
    va_end(v1);
}
*/

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
