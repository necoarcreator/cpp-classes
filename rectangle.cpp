
#include <iostream>
using namespace std;

class rectangle
{
    float length = 0, heigth = 0;
    float leftBelowPoint[2];
public:
    rectangle(float a, float b, float c, float d) : length{a}, heigth{b}, leftBelowPoint {c,d}
    {
    }

    float* getProperties()
    {
        float* M = new float[4];
        M[0] = length;
        M[1] = heigth;
        M[2] = leftBelowPoint[0];
        M[3] = leftBelowPoint[1];
        return M;
    }

    ~rectangle() {
        length = 0;
        heigth = 0;
        fill(begin(leftBelowPoint), end(leftBelowPoint), 0);
    };

};

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
int main()
{
    rectangle recA{ 2, 2, 1, 1};
    rectangle recB{ 1,1,5,5 };
    if (isCross(recA, recB))
    {
        cout << "Those rectangles are crossing" << endl;
    }
    else
    {
        cout << "Those rectangles aren't crossing" << endl;
    }
    }
    
