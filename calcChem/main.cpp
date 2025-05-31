#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <exception>
#include <sstream>
#include <algorithm>
#include <numeric>

using namespace std;
double molDole(vector<double> massWei, vector<double> masses, size_t i) //masses to moles
{
	double summ = 0.;
	for (size_t i = 0; i < massWei.size(); i++)
	{
		summ += massWei[i] / masses[i];
	}
	return massWei[i] / masses[i] / summ;
}
double massDole(vector<double> molWei, vector<double> masses, size_t i) //moles to masses
{
	double summ = 0.;
	for (size_t i = 0; i < molWei.size(); i++)
	{
		summ += molWei[i] * masses[i];
	}
	return molWei[i] * masses[i] / summ;
}
double molSummMass(vector<double> masses, vector<double> doles, bool isMolar)
{
	double summ = 0.;
	if (isMolar)
	{
		for (size_t i = 0; i < doles.size(); i++)
			{
			summ += doles[i] / masses[i];
			}
		return 1 / summ;
	}
	else
	{
		for (size_t i = 0; i < doles.size(); i++)
		{
			summ += doles[i] * masses[i];
		}
		return summ;
	}
}

bool findElem(string input, double* molWeightI)
{
	string line;
	map<string, double> molWeight;
	ifstream in("chem.txt"); // окрываем файл для чтения
	if (in.is_open())
	{
		while (getline(in, line))
		{
			stringstream ss(line);
			char delim = ' ';
			char delim2 = '/n';
			string x, y;
			getline(ss, x, delim);
			getline(ss, y, delim2);
			molWeight[x] = stod(y);
		}
	}
	else
	{
		cerr << "Can't open file to read";
		exit(1);
	}
	in.close();     // закрываем файл

	auto y = find_if(molWeight.begin(), molWeight.end(), [&](auto const& n)
		{
			return input == n.first;
		});
	if (y != molWeight.end())
	{
		*molWeightI = y->second;
		return true;
	}
	return false;
}

void input(size_t* n, vector<double>* output)
{
	string input;
	cout << "Insert number of elements" << endl;
	getline(cin, input);
	*n = static_cast<size_t>(stoi(input));

	bool isMolar;
	cout << "Choose if your case is molar dole (1) to usual or conversely (2)" << endl;
	getline(cin, input);
	if (input == "1" || input == "1 " || input ==  " 1") isMolar = true;
	else if (input == "2" || input ==  "2 " || input ==  " 2") isMolar = false;

	vector<double> doles;
	vector<double> masses;
	for (size_t i = 0; i < *n; i++)
	{
		cout << "Choose the name of elem according to its' name in Mendeleev's table" << endl;
		try
		{
			getline(cin, input);
			double* molVeightI = new double;
			if (findElem(input, molVeightI) == false) throw invalid_argument("No such elem in table");
			else masses.push_back(*molVeightI);
			delete molVeightI;
		}
		catch (const std::exception& e)
		{
			std::cerr << e.what() << std::endl;
			i--; 
		}
		if (isMolar) cout << "Choose the molar dole of elem" << endl;
		else cout << "Choose the mass dole of elem" << endl;
		getline(cin, input);
		doles.push_back(stod(input));

	}

	if (isMolar)
	{
		for (size_t j = 0; j < *n; j++)
		{
			output->push_back(massDole(doles, masses, j));
			cout << "w_i of elem " << j << " is " << (*output)[j] << endl;
		}

	}
	else
	{
		for (size_t j = 0; j < *n; j++)
		{
			output->push_back(massDole(doles, masses, j));
			cout << "x_i of elem " << j << " is " << (*output)[j] << endl;
		}	
	}
	cout << "Molar mass is " << molSummMass(masses, doles, isMolar);
	return;

}
void isoterm(bool isId, double T)
{
	size_t numAtm = 1000;
	double R = 8.31;

	double T_crit = 133;
	size_t type = 0;
	if (T == T_crit) type = 1;
	else if (T < T_crit) type = 2;
	
	vector<double> V;
	
	for (size_t i = 0; i < numAtm; i++)
	{
		V.push_back(i/1e5 + 4e-5);
	}
	vector<double> p = vector<double>(numAtm, 0.);

	string input;
	cout << "Insert number of elements" << endl;
	getline(cin, input);
	size_t n = static_cast<size_t>(stoi(input));
	
	vector<double>  weight = vector<double>(n);
	vector<double>  a = vector<double>(n);
	vector<double>  b = vector<double>(n);
	ifstream in("molesAir.txt"); // окрываем файл для чтения
	if (in.is_open())
	{
		size_t i = 0;
		while (getline(in, input) && i < n)
		{
			stringstream ss(input);
			char delim = ' ';
			char delim2 = '/n';
			string name, dole, aa, bb;
			getline(ss, name, delim);
			getline(ss, dole, delim);
			getline(ss, aa, delim);
			getline(ss, bb, delim2);
			weight[i] = stod(dole);
			a[i] = stod(aa);
			b[i] = stod(bb);
			i++;
		}
	}
	else
	{
		cerr << "Can't open file to read";
		exit(1);
	}
	in.close();     // закрываем файл
	if (!isId)
	{
		
		for (size_t i = 0; i < numAtm; i++)
		{
			double p_j = 0.;
			for (size_t j  = 0; j < n; j++)
			{
				p_j = R * T / (V[i] - b[j]) - a[j] / pow(V[i], 2);
				p[i] += p_j * weight[j]/100;
			}
		}

	}
	else
	{
		for (size_t i = 0; i < numAtm; i++)
		{
			double p_j = 0.;
			for (size_t j = 0; j < n; j++)
			{
				p_j = R * T / (V[i]);
				p[i] += p_j * weight[j]/100;
			}
		}
	}

	ofstream myfile;
	if (isId) myfile.open("p_id3.txt");
	else myfile.open("p_real22.txt");

	for (size_t i = 0; i < numAtm; i++)
	{
		myfile << V[i] << " " << p[i] << endl;
	}
	myfile.close();
}

int main(int argc, char* argv[])
{
	//size_t* N = new size_t;
	//vector<double>* weights = new vector<double>;
	//input(N, weights);
	

	bool isId;
	string input;
	cout << "Is your gas ideal? (1) - yes, (2) - no" << endl;
	getline(cin, input);
	if (input == "1" || input == "1 " || input == " 1") isId = true;
	else if (input == "2" || input == "2 " || input == " 2") isId = false;
	cout << "Choose temperature in K" << endl;
	getline(cin, input);

	isoterm(isId, stod(input));



	//delete N;
	//weights->clear();
	//delete weights;
}