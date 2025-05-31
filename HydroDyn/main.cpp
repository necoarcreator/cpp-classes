#include <iostream>
#include "godunov2.h"
#include <vector>

using namespace std; 

int main () {
	int test;
	cout << "Test # ";
	cin >> test;
	if (test > 6 || test < 0) {
		cout << "oops, that test doesn't exist" << endl;
		return 0;
	}
	/*
	int art_viscos = 0;
	cout << "Режим исскуственной вязкости" << endl;
	cout << "0 - без неё" << endl;
	cout << "1 - та, что была изначально" << endl;
	cout << "2 - константная" << endl;
	cout << "3 - линейная" << endl;
	cout << "4 - смешанная" << endl;
	cout << "5 - неймановская" << endl;
	cout << "6 - квадратичная" << endl;

	cin >> art_viscos;
	*/
	int isGodunov;
	cout << "Select a numeric method" << endl;
	cout << "1 - Godunov" << endl;
	cout << "0 - van Leer" << endl;
	cout << "2 - WENO" << endl;
	cout << "other - MacCormack" << endl;
	cin >> isGodunov;

	int min_mod_type=0;
	if (isGodunov == 0) {
		while (min_mod_type != 1 && min_mod_type != 2 && min_mod_type != 3 && min_mod_type != 4) {
			cout << "minmod type?" <<endl;
			cout << "1972 - 1" <<endl;
			cout << "1975 - 2" <<endl;
			cout << "1984 - 3" <<endl;
			cout << "'em all - 4" <<endl;
			cin >> min_mod_type;
		}
	}


	string belka;
	string belka2;
	string belka_1972;
	string belka_1975;
	string belka_1984;
	string belka_an;
	if (min_mod_type != 4) {
		cout << "Repository for results saving? ";
		cin >> belka;
		if (stoi(belka) > 9 || stoi(belka) < 0) {
			belka = "_new";
		}
		if (test == 0 || test == 1) {
			cout << "Repository for analytics saving? ";
			cin >> belka2;
		}
		else belka2 = "10";

		if (stoi(belka2) > 9 || stoi(belka2) < 0) {
			belka2 = "_new";
		}
	} else {
		cout << "Repository for results of 1972 saving? ";
		cin >> belka_1972;
		
		cout << "Repository for results of 1975 saving? ";
		cin >> belka_1975;
		
		cout << "Repository for results of 1984 saving? ";
		cin >> belka_1984;
		
		cout << "Repository for analytics saving? ";
		cin >> belka_an;
	}
	
	godunov god("input.txt");

	if (isGodunov == 1) god.Solver_godunov(test, "results" + belka + "/God", "results" + belka2 + "/Theory", ".csv", ',', 500);
	else if (isGodunov == 0) {
		if (min_mod_type == 1 || min_mod_type == 2 || min_mod_type == 3) {
		god.Solver_godunov_kolgan(test, "results" + belka + "/God", "results" + belka2 + "/Theory", ".csv", ',', 500, min_mod_type);
		} else {
			god.Solver_godunov_kolgan(test, "results" + belka_1972 + "/God", "results_new/Theory", ".csv", ',', 500, 1);
			god.Solver_godunov_kolgan(test, "results" + belka_1975 + "/God", "results_new/Theory", ".csv", ',', 500, 2);
			god.Solver_godunov_kolgan(test, "results" + belka_1984 + "/God", "results" + belka_an + "/Theory", ".csv", ',', 500, 3);
		}
	}
	else if (isGodunov == 2) god.Solver_godunov_WENO(test, "results" + belka + "/God", "results" + belka2 + "/Theory", ".csv", ',', 50);
	else god.Solver_mccormak(test, "results" + belka + "/God", "results" + belka2 + "/Theory", ".csv", ',', 50);


	return 0;
}
