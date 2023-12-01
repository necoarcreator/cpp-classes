#include "Queue.h"
#include "Stack.h"
#include "List.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <typeinfo> 

using namespace std;

template <typename T1>

T1 myReverse(T1& a)
{
    return a.rev(); //здесь я вшил разные методы в сами классы
}

template <typename T2>

void see(T2& a)
{
    a.print();
}

template <typename T3>

void reverseTZ(T3& a)
{
    if (a.isEmpty())
    {
        cout << "Error: the object is empty!\n";
        exit(1);
    }
    vector<int> Aint{ a.getA() };
    reverse(Aint.begin(), Aint.end());
    T3 b{ Aint };
    a = b;
    return;
}

template <>
void reverseTZ(List& a)
{
    if (a.isEmpty())
    {
        cout << "Error: the object is empty!\n";
        exit(1);

    }
    Node* current = a.head;
    Node* prev = NULL;
    Node* next = NULL;

    while (current != NULL) {
        next = current->next;
        current->next = prev;
        prev = current;
        current = next;
    }
    a.head = prev;
    cout << "The list has been reversed!\n";

}

int main()
{

    List A2;
    A2.push(2);
    A2.push(4);
    A2.push(8);
    A2.push(11);
    
    cout << "Reversing list:\n";

    see<List>(A2);

    reverseTZ<List>(A2);

    see<List>(A2);

    cout << "_______________________\n";
    cout << "Reversing queue:\n";

    Queue A({ 1, 2, 3 , 4 });

    reverseTZ<Queue>(A);

    cout << "_______________________\n";
    cout << "Reversing stack:\n";

    Stack B({ 11, 22, 33 , 44 });

    reverseTZ<Stack>(B);

}