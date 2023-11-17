#include <iostream>
#include <array>
#include <string>
#include <random>
#include <vector>
#include <math.h>
#include <limits.h>
#include <thread>
#include <chrono>
using namespace std;
int direct(int sideTr, int sideBull);
int finDir(int sideBV, int sideBH, int sideTV, int sideTH, int sideHV, int sideHH);
void fill(vector<int> &coordshB, vector<int> &coordsvB, vector<int>& coordshT, vector<int>& coordsvT, vector<int>& coordshH, vector<int>& coordsvH, array<array<string, 8>, 8> &pole);
void delk(vector<int>& coordshB, vector<int>& coordsvB, vector<int>& coordshT, vector<int>& coordsvT, vector<int>& coordshH, vector<int>& coordsvH, array<array<string, 8>, 8>& pole);
int random() {
    random_device rd; // non-deterministic generator
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_int_distribution<> dist(0, 7); // distribute results between 0 and 7 inclusive.
    return dist(gen); // pass the generator to the distribution.

}

int main()
{
    int m = 0, n = 0, k = 0;

    cout << "Enter the number of bulldozers: ";
    cin >> n;
    cout << "Enter the number of traaaashhhh: ";
    cin >> m;
    cout << "Enter the number of holes: ";
    cin >> k;

    int a = 8, b = 8;
    array<array<string, 8>, 8> pole;

    for (int i = 0; i < a; i++)
    {
        for (int j = 0; j < b; j++)
        {
            pole[i][j] = "|_|";
        }
    }

    vector<int> coordshB, coordshT, coordshH;
    vector<int> coordsvB, coordsvT, coordsvH;

    for (int i = 0; i < n; i++)
    {
        coordshB.push_back(0);
        coordsvB.push_back(0);
    }
    for (int i = 0; i < m; i++)
    {
        coordshT.push_back(0);
        coordsvT.push_back(0);
    }
    for (int i = 0; i < k; i++)
    {
        coordshH.push_back(0);
        coordsvH.push_back(0);
    }

    fill(coordshB, coordsvB, coordshT, coordsvT, coordshH, coordsvH, pole);

    for (int i = 0; i < a; i++)
    {
        for (int j = 0; j < b; j++)
        {
            cout << pole[i][j];
        }
        cout << endl;
    }
    cout << "----------------------------------" << endl;
    delk(coordshB, coordsvB, coordshT, coordsvT, coordshH, coordsvH, pole);
}

void fill(vector<int> & coordshB, vector<int> & coordsvB,
    vector<int>& coordshT, vector<int>& coordsvT,
    vector<int>& coordshH, vector<int>& coordsvH, array<array<string, 8>, 8> &pole)
{
    int n = coordshB.size();
    int m = coordshT.size();
    int k = coordshH.size();

    string obj = "|D|";
   
    for (int i = 0; i < n; i++)
    {
        coordshB[i] = random();
        coordsvB[i] = random();
        for (int j = 0; j < i; j++)
        {
            while ((coordshB[i] == coordshB[j]) and (coordsvB[i] == coordsvB[j]))
            {
                coordshB[j] = random();
                coordsvB[j] = random();
            }
        }
    }

    for (int l = 0; l < n; l++)
    {
        for (int i = 0; i < 8; i++)
        {

            for (int j = 0; j < 8; j++)
            {
                if ((i == coordshB[l]) and (j == coordsvB[l]))
                {
                    pole[i][j] = obj;
                }
            }
        }
    }
    obj = "|*|";

    for (int i = 0; i < m; i++)
    {
        coordshT[i] = random();
        coordsvT[i] = random();
        for (int j = 0; j < i; j++)
        {   
            if (i < n)
            {
                while (((coordshT[i] == coordshT[j]) and (coordsvT[i] == coordsvT[j])) || ((coordshB[i] == coordshT[j]) and (coordsvB[i] == coordsvT[j])))
                {
                    coordshT[j] = random();
                    coordsvT[j] = random();
                }
            }
            else
            {
                while ((coordshT[i] == coordshT[j]) and (coordsvT[i] == coordsvT[j]))
                {
                    coordshT[j] = random();
                    coordsvT[j] = random();
                }
            }
        }
    }
    for (int l = 0; l < m; l++)
    {
        for (int i = 0; i < 8; i++)
        {

            for (int j = 0; j < 8; j++)
            {
                if ((i == coordshT[l]) and (j == coordsvT[l]))
                {
                    pole[i][j] = obj;
                }
            }
        }
    }
    obj = "|0|";
    for (int i = 0; i < k; i++)
    {
        coordshH[i] = random();
        coordsvH[i] = random();
        for (int j = 0; j < i; j++)
        {   
            if ((i < n) and (i < m))
            {
                while (((coordshH[i] == coordshH[j]) and (coordsvH[i] == coordsvH[j])) || ((coordshB[i] == coordsvH[j]) and (coordsvB[i] == coordsvH[j])) || ((coordshT[i] == coordsvH[j]) and (coordsvT[i] == coordsvH[j])))
                {
                    coordshH[j] = random();
                    coordsvH[j] = random();
                }
            }
            else if (i < n)
            {
                while (((coordshH[i] == coordshH[j]) and (coordsvH[i] == coordsvH[j])) || ((coordshB[i] == coordsvH[j]) and (coordsvB[i] == coordsvH[j])))
                {
                    coordshH[j] = random();
                    coordsvH[j] = random();
                }
            }
            else if (i < m)
            {
                while (((coordshH[i] == coordshH[j]) and (coordsvH[i] == coordsvH[j])) || ((coordshT[i] == coordsvH[j]) and (coordsvT[i] == coordsvH[j])))
                {
                    coordshH[j] = random();
                    coordsvH[j] = random();
                }
            }
            else
            {
                while ((coordshH[i] == coordshH[j]) and (coordsvH[i] == coordsvH[j]))
                {
                    coordshH[j] = random();
                    coordsvH[j] = random();
                }
            }
        }
    }
    for (int l = 0; l < k; l++)
    {
        for (int i = 0; i < 8; i++)
        {

            for (int j = 0; j < 8; j++)
            {
                if ((i == coordshH[l]) and (j == coordsvH[l]))
                {
                    pole[i][j] = obj;
                }
            }
        }
    }
    return;


}

    int direct(int sideTr, int sideBull)
    {
        if (sideTr - sideBull > 1)
        {
            return 1; //бульдозер движ вверх/вправо
        }
        else if (sideTr - sideBull < 1)
        {
            return -1; //бульдозер движ вниз/влево
        }
        else
        {
            return 0; //не движ
        }
    }

    int finDir(int sideBV, int sideBH, int sideTV, int sideTH, int sideHV, int sideHH)
    {
        
        //нужна чтоб окончательно соориентировать бульдозер
        if ((sideHV - sideBV > 0) and (sideTV - sideBV > 0))
        {
            if (sideTH - sideBH > 0)
            {
                sideBH += 1;
            }
            else if (sideTH - sideBH < 0)
            {
                sideBH -= 1;
            }
        }
        return 0;
    }
//cout << coordshH.back() << "\n";
//cout << typeid(coordshH.back()).name();

void delk(vector<int>& coordshB, vector<int>& coordsvB,
    vector<int>& coordshT, vector<int>& coordsvT,
    vector<int>& coordshH, vector<int>& coordsvH, array<array<string, 8>, 8>& pole)
{
    int n = coordshB.size();
    int verT = coordsvT.back();
    int horT = coordshT.back();
    int k = coordshH.size();

    int lenmin = INT_MAX, ind = INT_MAX;

    for (int i = 0; i < n; i++)
    {
        if (abs(verT - coordsvB[i]) + abs(horT - coordshB[i]) < lenmin)
        {
            lenmin = abs(verT - coordsvB[i]) + abs(horT - coordshB[i]);
            ind = i;
        }
    } //находим ближайший дозер

    int hcount = 0;
    int dirV = direct(verT, coordsvB[ind]);
    int dirH = direct(horT, coordshB[ind]);
    bool isOnHole = false;
    //сближаемся
    while (dirV)
    {   
        if (isOnHole)
        {
            pole[coordshB[ind]][coordsvB[ind] - dirV] = "|0|";
            isOnHole = false;
        }

        if (pole[coordshB[ind]][coordsvB[ind]] == "|0|")
        {
            isOnHole = true;
        }
       pole[coordshB[ind]][coordsvB[ind]] = "|_|";
       coordsvB[ind] += dirV;
       hcount += 1;
       pole[coordshB[ind]][coordsvB[ind]] = "|D|";
       dirV = direct(verT, coordsvB[ind]);

       for (int i = 0; i < 8; i++)
       {
           for (int j = 0; j < 8; j++)
           {
               cout << pole[i][j];
           }
           cout << endl;
           
       }
        cout << "-------------------" << endl;
       this_thread::sleep_for(chrono::seconds(1));
    }
    while (dirH)
    {
        if (isOnHole)
        {
            pole[coordshB[ind] - dirH][coordsvB[ind]] = "|0|";
            isOnHole = false;
        }

        if (pole[coordshB[ind]][coordsvB[ind]] == "|0|")
        {
            isOnHole = true;
        }
        pole[coordshB[ind]][coordsvB[ind]] = "|_|";
        coordshB[ind] += dirH;
        hcount += 1;
        pole[coordshB[ind]][coordsvB[ind]] = "|D|";
        dirH = direct(horT, coordshB[ind]);

        for (int i = 0; i < 8; i++)
        {
            for (int j = 0; j < 8; j++)
            {
                cout << pole[i][j];
            }
            cout << endl;
            
        }
        cout << "-------------------" << endl;
        this_thread::sleep_for(chrono::seconds(1));
    }

    int lenmin2 = INT_MAX, indH = INT_MAX;

    for (int i = 0; i < k; i++)
    {
        if (abs(verT - coordsvH[i]) + abs(horT - coordshH[i]) < lenmin2)
        {
            lenmin2 = abs(verT - coordsvH[i]) + abs(horT - coordshH[i]);
            indH = i;
        }
    } //находим ближайшее о4о
    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            if ((j == coordsvH[indH]) and (i == coordshH[indH]))
            {
                cout << "|R|";
            }
            else
            {
                cout << pole[i][j];
            }
        }
        cout << endl;

    }
    //берем последнюю в списке кучу !
    //ищем ближайший к ней бульдозер !
    //сближаемся (поле рисуем) !
    //ищем ближайшее очо !
    //едем к нему
    //удаляем кучу
    coordsvT.pop_back();
    coordshT.pop_back();
    return;
}
