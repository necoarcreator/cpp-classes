
#include "Stack.h"
#include "Queue.h"
#include <iostream>

int main()
{
    Stack st({ 1,2,3 });

    st.push(4);

    cout << st.pop() << endl;



    Queue qe({ 5,6,7 });

    qe.push(8);

    cout << qe.pop();


}