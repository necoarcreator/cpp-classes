
#include <iostream>
#include <string>
#include <algorithm>
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
        if ((x[i] != '.') && (x[i] < '0' || x[i]>'9'))
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

unsigned long long int fact(int Nd)
{
	unsigned long long i = 1, sum = 1, Nl = Nd;

	if (Nl == 0)
	{
		return 1;
	}
	else
	{
		for( i = 1; i < Nl; i++)
		{
			sum *= i;
		}

		return sum;
	}
}


int ins(string N, bool isNorm)
{
    float Nf;
    int Nd;
    string Nnew;

    if (isNum(N))
    {
        cout << "Use numbers only" << std::endl; //по какой-то причине программа прекращает работу, если имеет дело со строчкой
        isNorm = true; //вернее, эту часть кода оона выполняет, и даже выводит на экран надпись о том, что пользователь - военный преступник, если вводит строку
    } // но при этом вылетает. Странно. Я ведб даже в отдельную переменную Nnew пишу всё в следующем вызове ins
    else
    {
        Nf = stof(N);
    }

    if ((not isNorm) && (Nf != round(Nf)))
    {
        cout << "Don't use fractional numbers" << std::endl; 
        isNorm = true;
    }
    else
    {
        Nd = stoi(N);
    }

    if ((not isNorm) && (Nd < 0))
    {
        cout << "Don't use negative numbers" << std::endl;
        isNorm = true;
    }

    if ((not isNorm) && (Nd == 0))
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

int main()
{
    int Nd;
	string Num;
	bool isNorm = false;
	std::cin >> Num;
	Nd = ins(Num, isNorm); //имею в вижу эту часть... Вот, тут Nnew. Конфликта записи данных быть,кажется, не должно. И запись результата ins в Nd невозможна, пока я не введу целое положительное...
	std::cout << fact(Nd) << std::endl;
	return 0;
}

