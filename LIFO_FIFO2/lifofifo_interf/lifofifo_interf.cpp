#include "Queue.h"
#include "Stack.h"
#include "List.h"
#include <iostream>
#include <vector>

using namespace std;

template <typename T1>

T1 reverse(T1 a)
{
    return a.rev(); //здесь я вшил разные методы в сами классы
}

template <typename T2>

void see(T2 a)
{
    a.print();
}


int main()
{
    
    //List<int> A({ 1, 2, 3 , 4});

    List<double> A2({ 1.4, 2.11, 3.33, 4.67});

    //List<long long int> A3({ 16368365856, 2257413616143, 34436457568678 , 43432424311 });

    //cout << A.search(3) << endl;

    //List<int> B(reverse(A));
    
    //see(A);
    //cout << B.search(3);
    //B.push(5);

    //B.pop();

}
