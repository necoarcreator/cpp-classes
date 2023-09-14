

#include <iostream>

int main()
{
    double val = 2;
    double* po = &val;
    double** p =  &po;
    std::cout << **p << std::endl;

    delete po;
    delete p;
    po = NULL;
    p = NULL;

    if (p == NULL)
    {
        std::cout << "Space cleared" << std::endl;
    }
    else
    {
        std::cout << p << std::endl;
        std::cout << **p << std::endl;

    }
}

