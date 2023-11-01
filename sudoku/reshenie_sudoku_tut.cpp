#include <vector>
#include <iostream>
#include <string>

using namespace std;

void helpToSort(int S, int n, vector<int> isk);
void findPairs(int S, int n, vector<int> all, vector<int> binV, int k);

int main()
{
	int S, n, size = 0, k;
	

	cout << "Enter the summ, number of elements in a group, forbidden numbers" << endl;
	cin >> S;
	cin >> n;
	vector<int> iskl;

	while (cin >> k && cin.peek() != '\n') //?
	{
		size += 1;
		iskl.push_back(k);
	}

	for (int j = 0; j < size; j++)
	{
		cout << "iskl[i] = " <<iskl[j] << endl;
	}

	helpToSort(S, n, iskl);

}

void findPairs(int S,int n, vector<int> all, vector<int> binV, int k)
{
	if (k == binV.size()) //at the end of the loop
	{
		int sum = 0;

		for (int l = 0; l < binV.size(); l++) // sum means how many numbers are included to array
		{
			sum += binV[l];
		}
			
		if (sum == n) //if it matches with the number of integers user needs
		{
			
			sum = 0; //now sum means the summ of numbers included

			for (int i = 0; i < all.size(); i++)
				sum += binV[i] * all[i];

			if (sum == S) //if it matches with the summ user needs
			{
				cout << "This combination may help you: ( | ";
				for (int j = 0; j < all.size(); j++)
				{	
					if ((binV[j]))
					{
						cout << all[j] << " | ";
					}
				}
				cout << ")" << endl;
				cout << endl;
				return;
			}
			else
			{
				return;
			}
		}
		else {
			return;
		}
		
	}
	
	
	
	
	binV[k] = 1; //try with k number first
	if (k != binV.size())
	{
		findPairs(S, n, all, binV, k + 1);
	}

	binV[k] = 0; //then try without it
	findPairs(S, n, all, binV, k + 1);

}

void helpToSort(int S, int n, vector<int> isk)
{
	vector <int> all{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };

	for (int i = 0; i < all.size(); i++)
	{
		for (auto k : isk)
		{
			if (all[i] == k)
			{
				all.erase(all.begin() + i);
			}
		}
	}

	vector <int> sootv;

	for (int j = 0; j < all.size(); j++)
	{
		sootv.push_back(0);
		
	}


	findPairs(S, n, all, sootv, 0);
	return;
}