#include "Queue.h"

Queue::Queue() : Stack() {};

Queue::Queue(vector<int> inpt) : Stack(inpt) {};

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