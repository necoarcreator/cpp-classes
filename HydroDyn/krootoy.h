#pragma once

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <fstream>
#include "tests.h"

using namespace std;

using Grid = vector <double>;

class krootoy {
	public:
	krootoy(string file_path);
	~krootoy();
	double Get_dt ();
	void makeGrid ();
	vector<string> pokushai (string file_path);
	
	void inletSlowTube(bool side); //дозвуковой в трубе
	
	void inletSlowEng(bool side);	//дозвуковой в двигателе
	
	void inletFast(bool side);	//сверхзвуковой
	
	void outletSlowPressure(bool side); //дозвуковой выход с заданием давления
	
	void freeFlow(bool side); //граница со свободным потоком
	
	void wallSlip(bool side); //стенка с проскальзыванием
	
	void wallSlipAdiab(bool side); //стенка с проскальзыванием и адиабатической температурой
	
	void wallNoSlip(bool side);	//стенка без проскальзывания
	
	void wallNoSlipLagrange();	//гу  для задачи с переменными Лагранжа

	void symmetry(bool side);	//условие симметрии
	
	void solverCrestNotDiv(string file_path, string format = ".csv", char _delim = '/t',
		size_t recordFreq = 5, int test = 1, int art_viscos = 0); //решалка)0

	void solverCrestDiv(string file_path, string format = ".csv", char _delim = '/t',
		size_t recordFreq = 5, int test = 1, int art_viscos = 0); //решалка в дивергентной форме
	
	void writeResults(string file_path, string format, char _delim, double _time); //я не буду писать это (вывод результатов)
	
	void ArtificialViscosity(); //искусственная вязкость - добавка к давлениямб определяемая градиентом скорости
	
	void ArtificialViscosityConst();

	void ArtificialViscosityLinear();

	void ArtificialViscosityNeumann();

	void ArtificialViscosityQuadro();

	void ArtificialViscosityMixed();

	double Lx_0;
	double Lx_1;
	double Ly_0;
	double Ly_1;
	double Lz_0;
	double Lz_1;
	double T_0;
	double T_1;
	double CFL;
	double C;
	size_t max_iter;
	size_t typeOfInlet, typeOfOutlet;
	
	double gamma;
	size_t fict;
	double h;
	int N_x, N_y, N_z;

	//in centers
	Grid* m; // Масса
	Grid* p; // Давление
	Grid* rho; //Плотность
	Grid* e; //Полная энергия
	Grid* ei; //Внутренняя энергия
	Grid* x; //Эйлерова координата
	Grid* v; //Скорость(только для иниализации)
	Grid* w; //Искусственная вязкость
	
	//in boundaries
	Grid* vb; //Скорость
	Grid* pb; //Давление
	Grid* xb; //Эйлерова координата

	//in centers
	Grid* p_new; //Давление
	Grid* rho_new; //Плотность
	Grid* e_new; //Полная энергия
	Grid* ei_new; // Внутренняя энергия
	Grid* m_new; // Масса
	
	//in boundaries
	Grid* vb_new; //Скорость
	Grid* xb_new; //Эйлерова координата

};

