#include <iostream>
using namespace std;

class A
{
    private:

    int a, b, c;

    public:

        

        A(int a1, int b1, int c1) : a(a1), b(b1), c(c1) {}

};

int main()
{
    A a1(1, 2, 3);
    A* a1_ptr = &a1;
    A** a1_ptr_ptr = &a1_ptr;
    int* k = (int*)a1_ptr;
    cout << "*a1_ptr = " << *a1_ptr << endl;
    cout << "*k = " << *k << endl;
    cout << "*(k + 1)= " << *(k + 1) << endl;
    cout << "*(k + 2)= " << *(k + 2) << endl;

    cout << "mem k = " << sizeof(k) << endl;
    cout << "mem (k + 1)= " << sizeof((k + 1)) << endl;
    cout << "mem (k + 2)= " << sizeof((k + 2)) << endl;

    cout << "mem a1_ptr = " << sizeof(a1_ptr) << endl;
    cout << "mem a1_ptr_ptr = " << sizeof(a1_ptr_ptr) << endl;
}