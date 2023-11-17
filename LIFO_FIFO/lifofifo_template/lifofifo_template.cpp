#include "Stack.h"
#include "Queue.h"


int main()
{
    vector <int> A{ 1,2,3,4,5};
    Stack inpt {-6, -10, -15, 15, 2};
    //Stack empt;
    inpt.push(6);
    int B = inpt.pop();

    cout << B;
    cout << endl;
    inpt.print();
    int C = inpt.pop();
    /*Queue inptt(vector<int>{1, 2, 3});
    inptt.push(6);
    int B = inptt.pop();

    cout << B;
    cout << endl;
    inptt.print();
    return 0;
    */
}
