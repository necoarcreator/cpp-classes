#include <vector>
#include <iostream>
#include <string>

using namespace std;

int helpPlayer(int S, int n, vector<int> isk);
unsigned long long int fact(int Nd);

int main()
{
	int S, n, size = 0, k;
	

	cout << "Enter the summ, number of elements in a group, forbidden numbers" << endl;
	cin >> S;
	cin >> n;
	vector<int> iskl;

	while (size < n && cin >> k && cin.peek() != '\n') //?
	{
		size += 1;
		iskl.push_back(k);
	}

	for (int j = 0; j < n; j++)
	{
		cout << "iskl[i] = " <<iskl[j] << endl;
	}

	helpPlayer(S, n, iskl);

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

string helpPlayer(int S, int n, vector<int> isk)
{
	vector <int> all{ 1, 2, 3, 4, 5, 6, 7, 8, 9 }; //?

	int len = pow(2, n) - 1; //musthave
	int j = 0;
	int sEff = S; //musthave
	int nEff = n; //musthave
	int* mid;
	string sum = "";

	while (j < isk.size()) //musthave
	{
		all.erase(isk[j]);
		sEff -= isk[j];
		nEff -= j;
		j += 1;
	};
	if (nEff == 1)
	{
		cout << sEff; //print lasted number
		return 0;
	}

	int combinations = fact(9 - j) / fact(nEff); //musthave

	float vzv = sEff / nEff; //?

	mid = &all[4];

	for (int i = 1; i < 9 - j; i++)
	{
		isk.push_back(i);
		sum += helpPlayer(Seff - i ; neff - 1; isk);

		cout << sum << endl;
	}

	

}
