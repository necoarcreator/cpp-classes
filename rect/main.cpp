
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "rectangle.cpp"

using namespace std;

bool isCross(rectangle A, rectangle B);

int main()
{
    rectangle recA{ 2, 2, 1, 1 };
    rectangle recB{ 1,1,5,5 };

    rectangle recC;

    recC.setProperties();
    if (isCross(recA, recC))
    {
        cout << "Those rectangles are crossing" << endl;
    }
    else
    {
        cout << "Those rectangles aren't crossing" << endl;
    }
}

bool isCross(rectangle A, rectangle B)
{
    float* aptr;
    float* bptr;

    aptr = A.getProperties();
    bptr = B.getProperties();

    cout << "aptr is " << endl;
    for (int i = 0; i < 4; i++)
    {
        cout << aptr[i] << "  ";
    }
    cout << endl;
    cout << "bptr is " << endl;
    for (int i = 0; i < 4; i++)
    {

        cout << bptr[i] << "  ";
    }
    cout << endl;

    if ((aptr[0] - bptr[0] - bptr[2]) > 0)
    {
        return false;
    }
    else if ((bptr[0] - aptr[0] - aptr[2]) > 0)
    {
        return false;
    }
    else if ((aptr[1] - bptr[1] - bptr[3]) > 0)
    {
        return false;
    }
    else if ((bptr[1] - aptr[1] - aptr[3]) > 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}
