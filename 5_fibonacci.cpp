
#include <iostream>
#include <vector>
#include <string>

using namespace std;

/*bool isNumeric(std::string const& str)
{
    return !str.empty() && std::all_of(str.begin(), str.end(), ::isdigit);
}
*/

bool isNum(string x)
{
    bool isN;
    int i = 0;

    if (x[0] == '-')
    {
        i += 1;
    }
    while (x[i])
    {
        if ((x[i] != '.')&&(x[i] < '0' || x[i]>'9'))
        { 
            isN = true; 
            break; 
        }
        i++;
    }
    if (isN)
    {
        return true;
    }
    else
    {
        return false;
    }
}


int ins(string N, bool isNorm)
{   
    float Nf;
    int Nd;
    string Nnew;

    if (isNum(N))
    {
        cout << "Use numbers only" << std::endl; //не знаю,почему новое считывание прекращается, хотя вот этот текст выдается
        isNorm = true;
    }
    else
    {
        Nf = stof(N);
    }

    if ((not isNorm)&&(Nf != round(Nf)))
    {
        cout << "Don't use fractional numbers" << std::endl; 
        isNorm = true;
    }
    else
    {
        Nd = stoi(N);
    }

    if ((not isNorm)&&(Nd < 0))
    {
        cout << "Don't use negative numbers" << std::endl;
        isNorm = true;
    }
      
    if ((not isNorm)&&(Nd == 0))
    {
        cout << "Don't use zero" << std::endl;
        isNorm = true;
    }
 

    if (isNorm)
    {
        std::cin >> Nnew;
        isNorm = false;
        Nd = ins(Nnew, isNorm);
    }

    return Nd;
}

int printNum(int N, vector<unsigned long long int> &fibs)
{
    int i;

    for (i = 0; i < N; i++)
    {
        cout << fibs[i] << " ";
    }
    return 0;
}

int getFibNum(int N)
{
    unsigned long long fib1 = 1, fib2 = 1, fib_sum = 0;
    vector<unsigned long long> fibs(N + 1);
    int i = 0;
    fibs[0] = 1;
    fibs[1] = 1;
    while (i < N - 2)
    {
        fib_sum = fib1 + fib2;
        fib1 = fib2;
        fib2 = fib_sum;
        fibs[i + 2] = fib_sum;
        i += 1;
    }
    printNum(N, fibs);
    return 0;
}

int main()
{
    string N;
    int Nd;
    bool isNorm = false;
    cin >> N;
    Nd = ins(N, isNorm);

    getFibNum(Nd);

}

