#include <iostream>
#include <ctime>

bool evenOdd(int n)
{
    if (n % 2)
    {
        return false;
    }
    else
    {
        return true;
    }
}

bool prove(int n)
{
    if (n % 3)
    {
        return false;
    }
    else
    {
        return true;
    }

}

bool evenOddBit(int n)
{
    int spr = 0b0001;

    if (spr & n)
    {
        return false;
    }
    else
    {
        return true;
    }


}

//#pragma optimize("", off)

int main()
{
    int i, j, n;
    long long N = 10e7;
    std::clock_t start;
    double duration_1, duration_2;

    srand(time(NULL));

    n = rand();
    start = std::clock();
    for (i = 0; i < N; i++)
    {
        n = rand();
        prove(n);
    }
    duration_1 = (std::clock() - start) / (double)CLOCKS_PER_SEC;

    n = rand();
    start = std::clock();
    for (i = 0; i < N; i++)
    {
        n = rand();
        evenOdd(n);
    }
    duration_2 = (std::clock() - start) / (double)CLOCKS_PER_SEC;

    std::cout << " The time for 2 is  " << duration_1 << std::endl;
    std::cout << " The time for 3 is  " << duration_2 << std::endl;

    std::cout << " The percent is  " << duration_2/duration_1 << "%" << std::endl;
}
//#pragma optimize("", on)