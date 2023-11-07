#include <iostream>
#include "Queue.h"
#include "Stack.h"
#include <vector>


using namespace std;

int main()
{
    Queue qeEmp;
    Queue qe({ 1,2,3 });

    Stack stEmp;
    Stack st({ 5,6,7 });

    qe.push(4);
    st.push(8);
    cout << qe.pop() << " is from queue\n";

    cout << st.pop() << " is from stack\n";

    cout << "\n let's try to pop something from a zero-sized queue:\n";

    qeEmp.pop();

    cout << "\n and now from a zero-sized stack:\n";

    stEmp.pop();
}
