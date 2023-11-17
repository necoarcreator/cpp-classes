
#include <iostream>
#include <array>
#include <string>
#include <random>/*
#include <time.h>
#include <thread>
#include <chrono>*/
#include <vector>

using namespace std;

void fill(int objType, int m, vector<int> &coordsh, vector<int> &coordsv, array<array<string, 8>, 8> &pole);

int random() {
    random_device rd; // non-deterministic generator
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_int_distribution<> dist(0, 7); // distribute results between 0 and 7 inclusive.

    /*
    srand(time(NULL));
    std::this_thread::sleep_for(std::chrono::seconds(1));
    int i = rand() % 8; */
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


    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            cout << pole[i][j];
        }
        cout << endl;
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

    fill(0, n, coordshB, coordsvB, pole);
    fill(1, m, coordshT, coordsvT, pole);
    fill(2, k, coordshH, coordsvH, pole);

    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            cout << pole[i][j];
        }
        cout << endl;
    }
}

void fill(int objType, int n, vector<int> &coordsh, vector<int> &coordsv, array<array<string, 8>, 8> &pole)
{
    string obj;
    int k = 0;
    if (objType == 0)
    {
        obj = "|D|";
    }
    else if (objType == 1)
    {
        obj = "|*|";
    }
    else
    {
        obj = "|0|";
    }

    for (int i = 0; i < n; i++)
    {
        coordsh[i] = random();
        coordsv[i] = random();
        for (int j = 0; j < i; j++)
        {
            while ((coordsh[i] == coordsh[j]) and (coordsv[i] == coordsv[j]))
            {
                coordsh[j] = random();
                coordsv[j] = random();
            }
        }
    }

    for (k = 0; k < n; k++)
    {
        for (int i = 0; i < 8; i++)
        {

            for (int j = 0; j < 8; j++)
            {
                if ((i == coordsh[k]) and (j == coordsv[k]))
                {
                    pole[i][j] = obj;
                }
            }
        }
    }

    return;


}
