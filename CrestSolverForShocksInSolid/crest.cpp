#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include "crest.h"

crest::crest(string file_path): krootoy(file_path) {
	s = new Grid(N_x, 0.0);
	s_new = new Grid(N_x, 0.0);
    epsilon = new Grid(N_x, 0.0);
    epsilon_new = new Grid(N_x, 0.0);

};

crest::~crest()  {
	delete s;
	delete s_new;
    delete epsilon;
    delete epsilon_new;
};

void crest::NewViscosity() {
	double mu_0 = 2.;
	for (ptrdiff_t i = fict; i < N_x - fict; i++) {
		if (((*vb)[i + 1] - (*vb)[i]) > 0) {
			(*w)[i] = (0.5 * mu_0 / (*rho)[i]) * ((*vb)[i + 1] - (*vb)[i]) * ((*vb)[i + 1] + (*vb)[i]);
		}

		else {
			(*w)[i] = 0;
		}
		(*w)[fict - 1] = (*w)[fict];
		(*w)[N_x - 1 - (fict - 1)] = (*w)[N_x - 1 - fict];
	}
	
	/*
	for (size_t i = 0; i < N_x; i++) {
		(*p)[i] = (*p)[i] + (*w)[i];
	}
	*/
	//cout << "Default" << endl;
}

void crest::Task1() {
	double p1, rho1, v1;
	double p2, rho2, v2;
	rho1 = 7850.0;
	for (ptrdiff_t i = 0; i < N_x; i++) {
		(*p)[i] = 0.0;
		(*rho)[i] = rho1;
		(*vb)[i] = 0.0;
		(*s)[i] = 0.0;
		(*epsilon)[i] = 0.0;
		(*v)[i] = 0.0;
	}
}

void crest::piston_no_wall(int v_0) {

	(*p)[0] = (*p)[1];
	(*rho)[0] = (*rho)[1];
	(*ei)[0] = (*ei)[1];
	//(*e)[fict - k] = (*e)[fict + k];
	(*m)[0] = (*m)[1];

	(*p)[N_x - fict] = (*p)[N_x - fict - 1];
	(*rho)[N_x - fict] = (*rho)[N_x - fict - 1];
	(*ei)[N_x - fict] = (*ei)[N_x - fict - 1];
		//(*e)[N_x - 1 - fict + k] = (*e)[N_x - fict - k - 1];
	(*m)[N_x - fict] = (*m)[N_x - fict - 1];

	(*vb)[0] = v_0;
	(*vb)[1] = v_0;

	(*vb)[N_x - fict] = (*vb)[N_x - fict - 1];
	
	/*
	for (ptrdiff_t k = 0; k <= fict; ++k) {
		(*p)[fict - k] = (*p)[fict + k];
		(*rho)[fict - k] = (*rho)[fict + k];
		(*ei)[fict - k] = (*ei)[fict + k];
		//(*e)[fict - k] = (*e)[fict + k];
		(*m)[fict - k] = (*m)[fict + k];

		(*p)[N_x - 1 - fict + k] = (*p)[N_x - fict - k - 1];
		(*rho)[N_x - 1 - fict + k] = (*rho)[N_x - fict - k - 1];
		(*ei)[N_x - 1 - fict + k] = (*ei)[N_x - fict - k - 1];
		//(*e)[N_x - 1 - fict + k] = (*e)[N_x - fict - k - 1];
		(*m)[N_x - 1 - fict + k] = (*m)[N_x - fict - k - 1];

		(*vb)[k] = v_0;
		(*vb)[N_x - 1 - k] = 0;
	}
	
	*/
	
	/*(*vb)[fict] = v_0;
	(*vb)[N_x - 1 - fict + 1] = 0;*/

}

void crest::piston_wall(int v_0) {

	(*p)[0] = (*p)[1];
	(*rho)[0] = (*rho)[1];
	(*ei)[0] = (*ei)[1];
	//(*e)[fict - k] = (*e)[fict + k];
	(*m)[0] = (*m)[1];

	(*p)[N_x - fict] = (*p)[N_x - fict - 1];
	(*rho)[N_x - fict] = (*rho)[N_x - fict - 1];
	(*ei)[N_x - fict] = (*ei)[N_x - fict - 1];
	//(*e)[N_x - 1 - fict + k] = (*e)[N_x - fict - k - 1];
	(*m)[N_x - fict] = (*m)[N_x - fict - 1];

	(*vb)[0] = v_0;
	(*vb)[1] = v_0;

	(*vb)[N_x - fict] = -(*vb)[N_x - fict - 1];

	/*
	for (ptrdiff_t k = 0; k <= fict; ++k) {
		(*p)[fict - k] = (*p)[fict + k];
		(*rho)[fict - k] = (*rho)[fict + k];
		(*ei)[fict - k] = (*ei)[fict + k];
		//(*e)[fict - k] = (*e)[fict + k];
		(*m)[fict - k] = (*m)[fict + k];

		(*p)[N_x - 1 - fict + k] = -(*p)[N_x - fict - k - 1];
		(*rho)[N_x - 1 - fict + k] = (*rho)[N_x - fict - k - 1];
		(*ei)[N_x - 1 - fict + k] = (*ei)[N_x - fict - k - 1];
		//(*e)[N_x - 1 - fict + k] = (*e)[N_x - fict - k - 1];
		(*m)[N_x - 1 - fict + k] = (*m)[N_x - fict - k - 1];

		(*vb)[k] = v_0;
		(*vb)[N_x - 1 - k] = -(*vb)[N_x - 1 - k - 1];
	}
	*/


	/*for (ptrdiff_t k = 1; k <= fict; ++k) {
		(*p)[fict - k] = (*p)[fict + k];
		(*rho)[fict - k] = (*rho)[fict + k];

		(*p)[N_x - 1 - fict + k] = -(*p)[N_x - fict - k - 1];
		(*rho)[N_x - 1 - fict + k] = -(*rho)[N_x - fict - k - 1];

		(*vb)[k] = v_0;
		(*vb)[N_x - 1 - k] = 0;
	}*/
}

void crest::solverCrestNotDiv(string file_path, string format, char _delim,
				size_t recordFreq, int test, int problem_type, int v_0) {
	
	makeGrid();
	
	switch(test) {
        case 1:
            Task1();
            break;
        default:
            std::cout << "Invalid test number" << std::endl;
            break;
    }

	for (ptrdiff_t i = fict; i < N_x - fict; i++) {
			(*ei)[i] = (*p)[i] / ((*rho)[i] * (gamma - 1));
			//(*e)[i] = (*ei)[i] + pow((*v)[i], 2) / 2;
			(*m)[i] = ((*xb)[i + 1] - (*xb)[i]) * (*rho)[i];
	}



	size_t k = 0;
	double dt = 0, timeAll = 0;

		while ((timeAll < (T_1 - T_0)) && (k < max_iter)) {
			
				switch(test) {
			case 1:
				piston_no_wall(v_0);
				break;
			case 2:
				piston_wall(v_0);
				break;
			default:
				std::cout << "Invalid test number" << std::endl;
				break;
		}
	
		//применяем вязкость
		NewViscosity ();
		
		dt = Get_dt();
		
		for (ptrdiff_t i = fict; i < N_x - fict + 1; i++) {
			// Считаем новое значение скорости на левой грани ячейки(поэтому массу берем среднюю)
			// По той же причине разностная прозводная p - левая
			(*vb_new)[i] = (*vb)[i] - dt * ((*p)[i] - (*p)[i - 1] + (*w)[i] - (*w)[i - 1]) / (0.5 * ((*m)[i] + (*m)[i - 1]));		
		}

		for (ptrdiff_t i = fict; i < N_x - fict + 1; i++) {
			// Считаем новое значение эйлеровы координаты
			(*xb_new)[i] = (*xb)[i] + dt * (*vb_new)[i];
		}
			
		for (ptrdiff_t i = fict; i < N_x - fict; i++) {
			// Считаем новое значение плотности координаты
			// Здесь мы берем правую разность, тк vb_new[i] - левая граница i - ой ячейки, vbnew[i + 1] - правая

			//(*rho_new)[i] = 1.0 / (1.0 / (*rho)[i] + dt * ((*vb_new)[i + 1] - (*vb_new)[i]) / (*m)[i]);
				
			// Так же очевидно, что тут можно посчитать плотность через новое значение эйлеровой координаты :
			// Формально это одно и то же, но, возможно, понятнее интуитивно
			// 
			(*rho_new)[i] = (*m)[i] / ((*xb_new)[i + 1] - (*xb_new)[i]);
		}
		
		for (ptrdiff_t i = fict; i < N_x - fict; i++) {
			// Считаем новое значение внутренней энергии - уравнение в НЕдивергентной форме
			//(*ei_new)[i] = (*ei)[i] / (1 + dt * (*rho_new)[i] * (gamma - 1) * ((*vb_new)[i + 1] - (*vb_new)[i]) / (*m)[i]);
			(*ei_new)[i] = (*ei)[i] + dt/(*m)[i] * ((*vb_new)[i]*((((*p)[i-1] * (*m)[i] + (*p)[i]*((*m)[i-1])) / ((*m)[i] + (*m)[i-1])) + 0.5 * ((*w)[i] + (*w)[i - 1])) - 

							((*vb_new)[i+1]*((((*p)[i+1] * (*m)[i] + (*p)[i]*((*m)[i+1])) / ((*m)[i] + (*m)[i+1])) + 0.5 * ((*w)[i] + (*w)[i + 1])))) +

							0.125 * (((*vb)[i+1] + (*vb)[i]) * ((*vb)[i+1] + (*vb)[i]) - ((*vb_new)[i+1] + (*vb_new)[i]) * ((*vb_new)[i+1] + (*vb_new)[i]));
		}
		
		for (ptrdiff_t i = fict; i < N_x - fict; i++) {
			// Считаем новое значение девиатора
			(*s_new)[i] = (*s)[i] + 2 * modul_sdwiga * (-(dt)*((*vb_new)[i+1] - (*vb_new)[i])/((*xb)[i+1]-(*xb)[i]) + 2.0/3.0 * (1.0/(*rho_new)[i] - 1.0/(*rho)[i]) / (1.0/(*rho_new)[i] + 1.0/(*rho)[i]));
		}
		for (ptrdiff_t i = fict; i < N_x - fict; i++) {
			// Считаем новое значение давления
			
			double p_hydr = modul_sjatia * (1 - rho_0/((*rho_new)[i]));
			
			if ((*s_new)[i] >= 2.0/3.0 * predel_teckuchesty) {
				(*p)[i] = p_hydr + 2.0 / 3.0 * (*s_new)[i]/abs((*s_new)[i]) * predel_teckuchesty;
			} else {
				(*p)[i] = p_hydr + (*s_new)[i];
			}
			
			(*p_new)[i] = (gamma - 1) * (*ei_new)[i] * (*rho_new)[i];
		}
		
		for (ptrdiff_t i = fict; i < N_x - fict; i++) {
			// Считаем новое значение деформации
			(*epsilon_new)[i] = (*epsilon)[i] + dt/(0.5 * ((*m)[i-1]) + (*m)[i]) * ((*vb)[i+1] - (*vb)[i-1]);
		}
		// Обновляем данные
		for (ptrdiff_t i = fict; i < N_x - fict; i++) {
			(*p)[i] = (*p_new)[i];
			(*rho)[i] = (*rho_new)[i];
			(*xb)[i] = (*xb_new)[i];
			(*vb)[i] = (*vb_new)[i];
			(*s)[i] = (*s_new)[i];
			(*epsilon)[i] = (*epsilon_new)[i];
			(*ei)[i] = (*ei_new)[i] + pow((((*vb)[i] + (*vb)[i + 1]) / 2),2) / 2;
		}
		if (k % recordFreq == 0) writeResults(file_path, format, _delim, timeAll);
		timeAll += dt;
		k++;
	}
	
	cout << "success" << endl;
	
	return;
}

void crest::writeResults (string file_path, string format = ".csv", char _delim = '/t', double _time = 0.) {

	for (ptrdiff_t i = fict; i < N_x - fict; i++) {
		(*x)[i] = ((*xb)[i] + (*xb)[i + 1]) / 2;
		(*v)[i] = ((*vb)[i] + (*vb)[i + 1]) / 2;
	}

	ofstream file;
	file.open (file_path + to_string(_time) + format);
	if (file.is_open()) {
		file << _time << endl;
		
		file << "x coordinate" << _delim << "x component of velocity" << _delim << "density" << _delim << "inner energy" << _delim << "pressure" << _delim << "deviator" << _delim << "deformation tens." << endl;
		
		for (ptrdiff_t i = 0; i < N_x; i++) {
			file << (*x)[i] << _delim << (*v)[i] << _delim << (*rho)[i] << _delim << (*e)[i] << _delim << (*p)[i] << _delim << (*s)[i] << _delim << (*epsilon)[i] << endl;
		}
	}
	else {
		cerr << "matvey durak" << endl;
		exit(1);
	}
	file.close();
}

/*
void krootoy::wallNoSlipLagrange() {
	(*p)[0] = (*p)[1];
	(*rho)[0] = (*rho)[1];
	(*p)[N_x - 1] = (*p)[N_x - 2];
	(*rho)[N_x - 1] = (*rho)[N_x - 2];
	(*vb)[1] = 0;
	(*vb)[N_x - 1] = 0;
}

void krootoy::symmetry(bool side = 0) {
	wallSlipAdiab(side);
}

double krootoy::Get_dt () {
	double dt_min = 1e6;
	double c_step = 0, dt = 0;	 

    for (int i = fict; i < N_x - fict; i++) {
		c_step = sqrt(gamma * (*p)[i] / (*rho)[i]);
		dt = ((*xb)[i + 1] - (*xb)[i]) / (abs((*vb)[i]) + c_step);

        if (dt < dt_min) {
			dt_min = dt;
		}
	}
    return dt_min * CFL;
}

void krootoy::makeGrid () {
	
	double dx_ = (Lx_1 - Lx_0) / (N_x - 2 * fict);
	 
	for (int i=0; i < N_x; i++) {
		(*x)[i] = (i + 0.5 - fict) * dx_;	
		(*xb)[i] = (i - fict) * dx_;
	}
	(*xb)[N_x] = (N_x - fict) * dx_; 
	
	return;
}

vector <string> krootoy::pokushai (string file_path) {
	ifstream file(file_path);
	vector <string> a;
	if(file.is_open())
	{
		string line, to_be_deleted;
		char delim1 = '\t', delim2 = '\n';
		int i = 0;
		while(getline(file, line)) {
			a.push_back("");
			stringstream sline(line);
			getline(sline, to_be_deleted, delim1);
			getline(sline, a[i], delim2);
			i++;
		}
	}
	else {
		cerr<< "A problem with opening file" << endl;
		exit(1);
	}
	file.close();
	return a;
}

void krootoy::ArtificialViscosity() {
	double mu_0 = 2.;
	for (size_t i = fict; i < N_x - fict; i++) {
		if (((*vb)[i + 1] - (*vb)[i]) < 0) {
			(*w)[i] = -mu_0 * (*rho)[i] * abs(((*vb)[i + 1] - (*vb)[i])) * ((*vb)[i + 1] - (*vb)[i]);
		}

		else {
			(*w)[i] = 0;
		}
		(*w)[0] = (*w)[1];
		(*w)[N_x - 1] = (*w)[N_x - 2];
	}

	for (size_t i = 0; i < N_x; i++) {
		(*p)[i] = (*p)[i] + (*w)[i];
	}
	//cout << "Default" << endl;
}
		
void krootoy::ArtificialViscosityConst() {
    double mu_0 = 2.;
    for (int i = 1; i < N_x - 1; ++i) {
        if ((*vb)[i + 1] - (*vb)[i] < 0.) {
            (*w)[i] = -mu_0 * (*m)[i];
        }
        else {
            (*w)[i] = 0.;
        }
    }
    (*w)[0] = (*w)[1];
    (*w)[N_x - 1] = (*w)[N_x - 2];
    for (int i = 0; i < N_x; ++i) {
        (*p)[i] = (*p)[i] + (*w)[i];
    }
	//cout << "Const" << endl;
}

void krootoy::ArtificialViscosityLinear() {
    double mu_0 = 2.;
    for (int i = 1; i < N_x - 1; ++i) {
        (*w)[i] = -mu_0 * (*rho)[i] * ((*vb)[i + 1] - (*vb)[i]);
    }
    (*w)[0] = (*w)[1];
    (*w)[N_x - 1] = (*w)[N_x - 2];
    for (int i = 0; i < N_x; ++i) {
        (*p)[i] = (*p)[i] + (*w)[i];
    }
	//cout << "Linear" << endl;
}

void krootoy::ArtificialViscosityNeumann() {
    double mu_0 = 2.;
    for (int i = 1; i < N_x - 1; ++i) {
        (*w)[i] = -mu_0 * (*rho)[i] * (*rho)[i] * abs((*vb)[i + 1] - (*vb)[i]) * ((*vb)[i + 1] - (*vb)[i]);
    }
    (*w)[0] = (*w)[1];
    (*w)[N_x - 1] = (*w)[N_x - 2];
    for (int i = 0; i < N_x; ++i) {
        (*p)[i] = (*p)[i] + (*w)[i];
    }
	//cout << "Neumann" << endl;
}

void krootoy::ArtificialViscosityQuadro() {
    double mu_0 = 2.;
    for (int i = 1; i < N_x - 1; ++i) {
        if ((*vb)[i + 1] - (*vb)[i] < 0.) {
            (*w)[i] = mu_0 * (*rho)[i] * ((*vb)[i + 1] - (*vb)[i]) * ((*vb)[i + 1] - (*vb)[i]);
        }
        else {
            (*w)[i] = 0.;
        }
    }
    (*w)[0] = (*w)[1];
    (*w)[N_x - 1] = (*w)[N_x - 2];
    for (int i = 0; i < N_x; ++i) {
        (*p)[i] = (*p)[i] + (*w)[i];
    }
	//cout << "quadro" << endl;
}

void krootoy::ArtificialViscosityMixed() {
    double mu_0 = 2.;
    for (int i = 1; i < N_x - 1; ++i) {
        if ((*vb)[i + 1] - (*vb)[i] < 0.) {
            (*w)[i] = mu_0 * (*rho)[i] * ((*vb)[i + 1] - (*vb)[i]) * ((*vb)[i + 1] - (*vb)[i]) - mu_0 * (*rho)[i] * ((*vb)[i + 1] - (*vb)[i]);
        }
        else {
            (*w)[i] = 0.;
        }
    }
    (*w)[0] = (*w)[1];
    (*w)[N_x - 1] = (*w)[N_x - 2];
    for (int i = 0; i < N_x; ++i) {
        (*p)[i] = (*p)[i] + (*w)[i];
    }
}

void krootoy::solverCrestNotDiv(string file_path, string format, char _delim,
				size_t recordFreq, int test, int art_viscos) {
	
	makeGrid();
	
	switch(test) {
        case 1:
            Test1(N_x, p, rho, v);
            break;
        case 2:
            Test2(N_x, p, rho, v);
            break;
        case 3:
            Test3(N_x, p, rho, v);
            break;
        case 4:
            Test4(N_x, p, rho, v);
            break;
        case 5:
            Test5(N_x, p, rho, v);
            break;
        default:
            std::cout << "Invalid test number" << std::endl;
            break;
    }

	for (size_t i = 0; i < N_x - fict; i++) {
			(*ei)[i] = (*p)[i] / ((*rho)[i] * (gamma - 1));
			(*e)[i] = (*ei)[i] + pow((*v)[i], 2) / 2;
			(*m)[i] = ((*xb)[i + 1] - (*xb)[i]) * (*rho)[i];
	}

	wallNoSlipLagrange();

	size_t k = 0;
	double dt = 0, timeAll = 0;

	while ((timeAll < (T_1 - T_0)) && (k < max_iter)) {
		
		switch(art_viscos) {
			case 1:	ArtificialViscosity(); break;
			case 2: ArtificialViscosityConst(); break;
			case 3: ArtificialViscosityLinear(); break;
			case 4: ArtificialViscosityMixed(); break;
			case 5: ArtificialViscosityNeumann(); break;
			case 6: ArtificialViscosityQuadro(); break;
			default: break;
		}
		
		dt = Get_dt();
		
		for (size_t i = fict; i < N_x - fict + 1; i++) {
			// Считаем новое значение скорости на левой грани ячейки(поэтому массу берем среднюю)
			// По той же причине разностная прозводная p - левая
			(*vb_new)[i] = (*vb)[i] - dt * ((*p)[i] - (*p)[i - 1]) / (.5 * ((*m)[i] + (*m)[i - 1]));		
		}

		for (size_t i = fict; i < N_x - fict + 1; i++) {
			// Считаем новое значение эйлеровы координаты
			(*xb_new)[i] = (*xb)[i] + dt * (*vb_new)[i];
		}
			
		for (size_t i = fict; i < N_x - fict; i++) {
			// Считаем новое значение плотности координаты
			// Здесь мы берем правую разность, тк vb_new[i] - левая граница i - ой ячейки, vbnew[i + 1] - правая

			(*rho_new)[i] = 1.0 / (1.0 / (*rho)[i] + dt * ((*vb_new)[i + 1] - (*vb_new)[i]) / (*m)[i]);
				
			// Так же очевидно, что тут можно посчитать плотность через новое значение эйлеровой координаты :
			// Формально это одно и то же, но, возможно, понятнее интуитивно
			// 
			// rho_new[i] = m[i] / (xb_new[i + 1] - xb_new[i])
		}
		
		for (size_t i = fict; i < N_x - fict; i++) {
			// Считаем новое значение внутренней энергии - уравнение в НЕдивергентной форме
			(*ei_new)[i] = (*ei)[i] / (1 + dt * (*rho_new)[i] * (gamma - 1) * ((*vb_new)[i + 1] - (*vb_new)[i]) / (*m)[i]);
		}

		for (size_t i = fict; i < N_x - fict; i++) {
			// Считаем новое значение давления
			(*p_new)[i] = (gamma - 1) * (*ei_new)[i] * (*rho_new)[i];
		}
					
		// Обновляем данные
		for (size_t i = fict; i < N_x - fict; i++) {
			(*p)[i] = (*p_new)[i];
			(*rho)[i] = (*rho_new)[i];
			(*xb)[i] = (*xb_new)[i];
			(*vb)[i] = (*vb_new)[i];
			(*ei)[i] = (*ei_new)[i];
			(*e)[i] = (*ei)[i] + pow((((*vb)[i] + (*vb)[i + 1]) / 2),2) / 2;
		}
		if (k % recordFreq == 0) writeResults(file_path, format, _delim, timeAll);
		timeAll += dt;
		k++;
	}
	
	cout << "success" << endl;
	
	return;
}

void krootoy::solverCrestDiv(string file_path, string format, char _delim,
					size_t recordFreq, int test, int art_viscos) {

	makeGrid();

	switch(test) {
        case 1:
            Test1(N_x, p, rho, v);
            break;
        case 2:
            Test2(N_x, p, rho, v);
            break;
        case 3:
            Test3(N_x, p, rho, v);
            break;
        case 4:
            Test4(N_x, p, rho, v);
            break;
        case 5:
            Test5(N_x, p, rho, v);
            break;
        default:
            std::cout << "Invalid test number" << std::endl;
            break;
    }

	for (size_t i = 0; i < N_x - fict; i++) {
		(*ei)[i] = (*p)[i] / ((*rho)[i] * (gamma - 1));
		(*e)[i] = (*ei)[i] + pow((*v)[i], 2) / 2;
		(*m)[i] = ((*xb)[i + 1] - (*xb)[i]) * (*rho)[i];
	}

	wallNoSlipLagrange();

	size_t k = 0;
	double dt = 0, timeAll = 0;

	while ((timeAll < (T_1 - T_0)) && (k < max_iter)) {
		
		switch(art_viscos) {
			case 1:	ArtificialViscosity(); break;
			case 2: ArtificialViscosityConst(); break;
			case 3: ArtificialViscosityLinear(); break;
			case 4: ArtificialViscosityMixed(); break;
			case 5: ArtificialViscosityNeumann(); break;
			case 6: ArtificialViscosityQuadro(); break;
			default: break;
		}
		
		dt = Get_dt();

		for (size_t i = fict; i < N_x - fict + 1; i++) {
			// с помощью линейной интерполяции нашли давление на левой грани ячейки
			(*pb)[i] = ((*p)[i] * (*m)[i - 1] + (*p)[i - 1] * (*m)[i]) / ((*m)[i] + (*m)[i - 1]);
			(*pb)[0] = (*pb)[1];
			(*pb)[N_x - 1] = (*pb)[N_x - 2];
		}
			
		for (size_t i = fict; i < N_x - fict + 1; i++) {
			// Считаем новое значение скорости на левой грани ячейки(поэтому массу берем среднюю)
			// По той же причине разностная прозводная p - левая

			(*vb_new)[i] = (*vb)[i] - dt * ((*p)[i] - (*p)[i - 1]) / (.5 * ((*m)[i] + (*m)[i - 1]));
		}
				
		for (size_t i = fict; i < N_x - fict + 1; i++) {
			// Считаем новое значение эйлеровы координаты
			(*xb_new)[i] = (*xb)[i] + dt * (*vb_new)[i];
		}
					
		for (size_t i = fict; i < N_x - fict; i++) {
			// Считаем новое значение плотности координаты
			// 
			// Здесь мы берем правую разность, тк vb_new[i] - левая граница i - ой ячейки, vbnew[i + 1] - правая
			// rho_new[i] = 1.0 / (1.0 / (*rho)[i] + dt * (vb_new[i + 1] - vb_new[i]) / m[i])
			//
			// Так же очевидно, что тут можно посчитать плотность через новое значение эйлеровой координаты :
			// Формально это одно и то же, но, возможно, понятнее интуитивно

			(*rho_new)[i] = (*m)[i] / ((*xb_new)[i + 1] - (*xb_new)[i]);
		}
			
		for (size_t i = fict; i < N_x - fict; i++) {
			// Считаем новое значение полной энергии - уравнение в дивергентной форме
			(*e_new)[i] = (*e)[i] - dt * ((*pb)[i + 1] * (*vb_new)[i + 1] - (*pb)[i] * (*vb_new)[i]) / (*m)[i];
		}
				
		
			for (size_t i = fict; i < N_x - fict; i++) {
				// Считаем новое значение внутренней энергии через полную энергию и скорость
				(*ei_new)[i] = (*e_new)[i] - pow((((*vb)[i] + (*vb)[i + 1]) / 2), 2) / 2;
			}
					
		for (size_t i = fict; i < N_x - fict; i++) {
			// Считаем новое значение давления
			(*p_new)[i] = (gamma - 1) * (*ei_new)[i] * (*rho_new)[i];
		}
						
		// Обновляем данные
		for (size_t i = fict; i < N_x - fict; i++) {
			(*p)[i] = (*p_new)[i];
			(*rho)[i] = (*rho_new)[i];
			(*xb)[i] = (*xb_new)[i];
			(*vb)[i] = (*vb_new)[i];
			(*ei)[i] = (*ei_new)[i];
			(*e)[i] = (*e_new)[i];
		}
		
		if (k % recordFreq == 0) writeResults(file_path, format, _delim, timeAll);
		timeAll += dt;
		k++;
	}

	cout << "success" << endl;

	
	return;
}


void krootoy::writeResults(string file_path, string format = ".csv", char _delim = '/t', double _time = 0.) {

	for (size_t i = fict; i < N_x - fict; i++) {
		(*x)[i] = ((*xb)[i] + (*xb)[i + 1]) / 2;
		(*v)[i] = ((*vb)[i] + (*vb)[i + 1]) / 2;
	}

	ofstream file;
	file.open (file_path + to_string(_time) + format);
	if (file.is_open()) {
		file << _time << endl;
		
		file << "x coordinate" << _delim << "x component of velocity" << _delim << "density" << _delim << "inner energy" << _delim << "pressure" << endl;
		
		for (size_t i = 0; i < N_x; i++) {
			file << (*x)[i] << _delim << (*v)[i] << _delim << (*rho)[i] << _delim << (*e)[i] << _delim << (*p)[i] << endl;
		}
	}
	else {
		cerr << "matvey durak" << endl;
		exit(1);
	}
	file.close();
}




// Toro tests

void Test1(size_t N_x, vector<double>* p, vector<double>* rho, vector<double>* v) {
	double p1, rho1, v1;
	double p2, rho2, v2;
	rho1 = 1.0;
	v1 = 0.0;
	p1 = 1.0;
	rho2 = 0.125;
	v2 = 0.0;
	p2 = 0.1;
	for (int i = 0; i < N_x; i++) {
		if (i < 0.5 * N_x) {
			(*p)[i] = p1;
			(*rho)[i] = rho1;
			(*v)[i] = v1;
		}
		else {
			(*p)[i] = p2;
			(*rho)[i] = rho2;
			(*v)[i] = v2;
		}
	}
}

void Test2(size_t N_x, vector<double>* p, vector<double>* rho, vector<double>* v) {
	double p1, rho1, v1;
	double p2, rho2, v2;
	rho1 = 1.0;
	v1 = -2.0;
	p1 = 0.4;
	rho2 = 1.0;
	v2 = 2.0;
	p2 = 0.4;
	for (int i = 0; i < N_x; i++) {
		if (i < 0.5 * N_x) {
			(*p)[i] = p1;
			(*rho)[i] = rho1;
			(*v)[i] = v1;
		}
		else {
			(*p)[i] = p2;
			(*rho)[i] = rho2;
			(*v)[i] = v2;
		}
	}
}

void Test3(size_t N_x, vector<double>* p, vector<double>* rho, vector<double>* v) {
	double p1, rho1, v1;
	double p2, rho2, v2;
	rho1 = 1.0;
	v1 = 0.0;
	p1 = 1000.0;
	rho2 = 1.0;
	v2 = 0.0;
	p2 = 0.01;
	for (int i = 0; i < N_x; i++) {
		if (i < 0.5 * N_x) {
			(*p)[i] = p1;
			(*rho)[i] = rho1;
			(*v)[i] = v1;
		}
		else {
			(*p)[i] = p2;
			(*rho)[i] = rho2;
			(*v)[i] = v2;
		}
	}
}

void Test4(size_t N_x, vector<double>* p, vector<double>* rho, vector<double>* v) {
	double p1, rho1, v1;
	double p2, rho2, v2;
	rho1 = 1.0;
	v1 = 0.0;
	p1 = 0.01;
	rho2 = 1.0;
	v2 = 0.0;
	p2 = 100.0;
	for (int i = 0; i < N_x; i++) {
		if (i < 0.5 * N_x) {
			(*p)[i] = p1;
			(*rho)[i] = rho1;
			(*v)[i] = v1;
		}
		else {
			(*p)[i] = p2;
			(*rho)[i] = rho2;
			(*v)[i] = v2;
		}
	}
}

void Test5(size_t N_x, vector<double>* p, vector<double>* rho, vector<double>* v) {
	double p1, rho1, v1;
	double p2, rho2, v2;
	rho1 = 5.99924;
	v1 = 19.5975;
	p1 = 460.894;
	rho2 = 5.99242;
	v2 = -6.19633;
	p2 = 46.0950;
	for (int i = 0; i < N_x; i++) {
		if (i < 0.5 * N_x) {
			(*p)[i] = p1;
			(*rho)[i] = rho1;
			(*v)[i] = v1;
		}
		else {
			(*p)[i] = p2;
			(*rho)[i] = rho2;
			(*v)[i] = v2;
		}
	}
}

*/