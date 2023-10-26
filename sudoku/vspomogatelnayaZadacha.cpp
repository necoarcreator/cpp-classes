#include <vector>
#include <iostream>
#include <string>

using namespace std;
vector<string> filt(vector <string> arr, string a);
vector<string> help(vector <string> arr, int len);
unsigned long long int fact(int Nd);

//vector<vector<string>> Reshape(vector<string> res, int numElem);
//vector<string> soch(vector<vector<string>> per, int numElem);

int main()
{
    vector<string> arr = { "a", "b", "c"};
    vector <string> res = help(arr, arr.size());
    int len = arr.size(), lenr = res.size();
    cout << endl;
    cout << len << " and" << lenr;
    for (int j = 0; j < lenr; j++)
    {
        cout << res[j] << " ";
        if (!((j + 1) % len))
        {
            cout << endl;
        }
        
    }
    /*
    vector <string> res = filt(arr, "c");
    int lenr = res.size();
    for (int j = 0; j < lenr; j++)
    {
        cout << res[j] << " ";
    }
    */
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
        for (i = 1; i < Nl; i++)
        {
            sum *= i;
        }

        return sum;
    }
}


vector<string> filt(vector <string> arr, string a)
{
    vector<string> res;
    for (string i : arr)
    {
        if (a.compare(i))
        {
            res.push_back(i);
        }
    }
    return res;
}

vector<string> help(vector <string> arr, int len)
{
    vector<string> result;

    if (arr.size() == 1)
    {
        return arr;
    }
    else {
        int numbvst = fact(len - 1);
        for (string let : arr)
        {
            //cout << "Now we check for" << let << endl;

            vector<string> temp = help(filt(arr, let), len);
           
            if (arr.size() > 1) {
                for (int kal = 0; kal < numbvst; kal++) {
                    temp.insert(temp.begin() + kal * arr.size(), let);
                }
            }
            
            result.insert(result.end(), temp.begin(), temp.end());
           
            
        }
    }
    return result;
}
/*
vector<vector<string>> Reshape(vector<string> res, int numElem) // тут я превращаю огромный вектор на входе в вектор векторов, его элементы - это перестановки размерности n
{
    vector<vector<string>> perest;
    int len = res.size();
    int lenper = len / numElem;
    for (int k = 0; i < lenper; i++)
    {
        for (int j = 0; j < numElem; j++)
            {
                perest[k][j] = res[k + j];
            }
    }

    return perest;
}

vector<string> soch(vector<vector<string>> per, int lenper, int numElem) //тут будет какой-то из алгоритмов сортировки внутри каждой из перестановок. Быстрая, пузырьком - неважно

{
    int numrow = lenper / numElem;

    for (int i = 0; i < numrow; i++) //по всем перест
    {
        for (int j = 0; j < )
    }


    //а после сортировки я просто сравниваю бывшие перестановки и совпадающие удаляю. Вот так и получем сочетания.



}
*/