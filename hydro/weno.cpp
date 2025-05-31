#include <iostream>
#include <fstream>
#include <string> 
#include <vector>
#include <map>
#include <stdlib.h>
#include <sstream>
#include <mpi.h>
#include <cassert>
#include <algorithm>
#include <numeric>
#define _USE_MATH_DEFINES 
#include <math.h>
typedef std::vector<double> vec;

template <typename T>
class Field2D {
public:
  Field2D(int nx, int ny) : nx_(nx), ny_(ny), data_((nx)* (ny)) {}

  T& operator()(int i, int j) {
    assert(i >= 0 && i < nx_);
    assert(j >= 0 && j < ny_);
    return data_[i + nx_ * j];
  }

  const T& operator()(int i, int j) const {
    assert(i >= 0 && i < nx_);
    assert(j >= 0 && j < ny_);
    return data_[i + nx_ * j];
  }

  T* data() { return data_.data(); }
  const T* data() const { return data_.data(); }

  int nx() const { return nx_; }
  int ny() const { return ny_; }
private:
  int nx_, ny_;
  std::vector<T> data_;
};

struct Params
{
	double CFL;					// число Куранта
	double gamma;				// показатель адиабаты
	double dx, dy;				// шаг по пространству
	int scheme_type;			// тип расчетной схемы: Годунов или Годунов-Колган или Годунов-Колган-Родионов
	double Lx, Ly;				// размер расчетной области
	double T;					// период расчета
	int Nx, Ny;					// число шагов
	double Quser;				// максимальный скачок давлений
	int fict;					// число фиктивных шагов

	Params()
	{
		this->CFL = 0.05;
		this->gamma = 1.4;
		this->Nx = 20;
		this->Lx = 1.0;
		this->Ly = 1.0;
		this->Ny = 20;
		this->T = 1.;
		this->Quser = 2.;
		this->dx = Lx / Nx;
		this->dy = Ly / Ny;
		this->scheme_type = 1;
		if (scheme_type == 0) this->fict = 1;
		else this->fict = 2;
	}

	Params(int Nx_, int Ny_, int type, double cfl_, double Lx_, double Ly_, double T_, double gamma_, double Quser_)
	{

		this->Nx = Nx_;// +2 * this->fict;
		this->Ny = Ny_;// +2 * this->fict;
		
		this->scheme_type = type;
		if (type == 0) this->fict = 1;
		else if (type == 1) fict = 2;
		else this->fict = 5;
		
		this->CFL = cfl_;
		
		this->Lx = Lx_;
		this->Ly = Ly_;
		
		this->T = T_;
	
		this->gamma = gamma_;
		
		this->Quser = Quser_;
		this->dx = Lx / Nx_;
		this->dy = Ly / Ny_;
	}
};

std::vector<int> delimetrs(int p) {
	std::vector<int> res(2, 0);
	std::vector<int> delimeters;

	if (p == 1) {
		res[0] = 1;
		res[1] = 1;
		return res;
	}

	// Собираем все делители
	for (int i = 1; i <= sqrt(p); ++i) {
		if (p % i == 0) {
			delimeters.push_back(i);
			if (i != p / i) {
				delimeters.push_back(p / i);
			}
		}
	}
	//if (p != 2 && delimeters.size() <= 2) {
	//	return delimetrs(p - 1);
	//}
	// Сортируем делители
	std::sort(delimeters.begin(), delimeters.end());

	int product1 = 1, product2 = p;
	int min_diff = p - 1; // Начальная минимальная разница

	// Ищем пару делителей с минимальной разницей
	for (int del : delimeters) {
		int a = del;
		int b = p / a;
		if (a * b == p) {
			int diff = abs(a - b);
			if (diff < min_diff) {
				min_diff = diff;
				product1 = a;
				product2 = b;
			}
		}
	}

	res[1] = product1;
	res[0] = product2;

	return res;
}


void log_ghost_cells(int rank, int step, Field2D<double> &rho, int local_Nx, int local_Ny, int fict) {
    std::cout << "Rank " << rank << " Step " << step 
              << " Left ghost rho(0,5): " << rho(0, 5) 
              << " Right ghost rho(" << local_Nx + fict << ",5): " << rho(local_Nx + fict, 5) << std::endl;
    std::cout << "Rank " << rank << " Step " << step 
              << " Left calc rho(" << fict << ",5): " << rho(fict, 5) 
              << " Right calc rho(" << local_Nx + fict - 1 << ",5): " << rho(local_Nx + fict - 1, 5) << std::endl;
}


void data_to_file(int local_Nx, int local_Ny, int fict, double gamma, double time, vec x, vec y, Field2D<double> p, 
					Field2D<double> vx, Field2D<double> vy, Field2D<double> rho, int rank, std::vector<int> dims, const char* name) {
	std::ostringstream bufr;
	for (int i = fict; i < local_Nx + fict; ++i) {
		for (int j = fict; j < local_Ny + fict; ++j) {

			bufr << x[i] << "," << y[j] << "," << p(i, j) << ","
				<< vx(i, j) << "," << vy(i, j) << "," << rho(i, j) << "," << p(i, j) / ((gamma - 1.0) * rho(i, j)) << "\n";
		}
	}

	std::string local_data = bufr.str();
	int local_size = local_data.size();

	int offset = 0;
	MPI_Exscan(
		&local_size, // Вход: размер данных текущего процесса
		&offset, // Выход: сумма размеров предыдущих процессов
		1, // Количество элементов
		MPI_INT, // Тип данных
		MPI_SUM, // Операция суммирования
		MPI_COMM_WORLD // Коммуникатор
	);

	char name_file[100];
	
	sprintf(name_file, "%sout_%f_.csv", name, time); // для кластера

	//sprintf_s(name_file, "%sout_%f_.csv", name, time); // для локального компа

	MPI_File fh;
	MPI_File_open(
		MPI_COMM_WORLD, // Все процессы
		name_file, // Имя файла
		MPI_MODE_CREATE | // Создать файл
		MPI_MODE_WRONLY, // Только запись
		MPI_INFO_NULL, // Доп. параметры
		&fh // Хэндлер файла
	);

	const char* header = "x,y,p,vx,vy,r,e;\n";
	int header_len = strlen(header);
	if (rank == 0) {
		MPI_File_write(fh, header, header_len, MPI_CHAR, MPI_STATUS_IGNORE);
	}
	
	offset += header_len;

	MPI_File_write_at_all(
		fh, // Хэндлер
		offset, // Смещение
		local_data.c_str(), // Данные
		local_size, // Размер
		MPI_CHAR, // Тип данных
		MPI_STATUS_IGNORE // Игнорировать статус
	);

	MPI_File_close(&fh);

}

// from conservative
void convert_from_conservative(double gamma, double& p, double& vx, double& vy, double& rho, double& m, double& impx, double& impy, double& e)
{
	if (m <= 0.0 || std::isnan(m)) {
		std::cerr << "‼️ [convert] rho = " << rho << ", m = " << m << std::endl;
	}

	p = (gamma - 1.0) * (e - 0.5 * (pow(impx, 2.0) + pow(impy, 2.0)) / m);
	vx = impx / m;
	vy = impy / m;
	rho = m;

	
	if (std::isnan(impx) || std::isnan(impy) || std::isnan(e)) {
		std::cerr << "‼️ [convert] NaN detected: "
			<< "impx = " << impx << ", impy = " << impy
			<< ", e = " << e << std::endl;
	}
}
// to conservative
void convert_to_conservative(double gamma, double& p, double& vx, double& vy, double& rho, double& m, double& impx, double& impy, double& e)
{
	if (rho <= 0.0 || std::isnan(rho)) {
		std::cerr << "‼️ [convert] rho = " << rho << ", m = " << m << std::endl;
	}

	m = rho;

	impx = rho * vx;
	impy = rho * vy;
	e = 0.5 * rho * (pow(vx, 2.0) + pow(vy, 2.0)) + p / (gamma - 1.0);

	if (std::isnan(p) || std::isnan(vx) || std::isnan(vy)) {
		std::cerr << "‼️ [convert] NaN detected: "
			<< "p = " << p << ", vx = " << vx << ", vy = " << vy
			<< ", e = " << e << ", rho = " << rho << std::endl;
	}
}


void boundary_cond_x(double gamma, double m, double impx, double impy, double e, double& mb, double& impxb, double& impyb, double& eb, int mode)
{
	double mach, rho1, p1, vx1, vy1, m1, impx1, impy1, e1;

	// wall
	if (mode == 0) {
		mb = m;
		impxb = -impx;
		impyb = impy;
		eb = e;
	}

	// free flux
	else if (mode == 1) {
		mb = m;
		impxb = impx;
		impyb = impy;
		eb = e;
	}

	//mach params
	else {
		mach = 3.0;
		rho1 = 1.225;

		p1 = 50.0 * 1e3;
		vx1 = mach * std::sqrt(gamma * p1 / rho1);
		vy1 = impy / rho1;
		convert_to_conservative(gamma, p1, vx1, vy1, rho1, m1, impx1, impy1, e1);

		mb = m1;
		impxb = impx1;
		impyb = 0.0;
		eb = e1;
	}

}

void boundary_cond_y(double gamma, double m, double impx, double impy, double e, double& mb, double& impxb, double& impyb, double& eb, int mode)
{
	double mach, rho1, p1, vx1, vy1, m1, impx1, impy1, e1;

	// wall
	if (mode == 0) {
		mb = m;
		impxb = impx;
		impyb = -impy;
		eb = e;
	}

	// free flux
	else if (mode == 1) {
		mb = m;
		impxb = impx;
		impyb = impy;
		eb = e;
	}

	//mach params
	else {
		mach = 3.0;
		rho1 = 1.225;

		p1 = 50.0 * 1e3;
		vx1 = impx / rho1;
		vy1 = mach * std::sqrt(gamma * p1 / rho1);
		convert_to_conservative(gamma, p1, vx1, vy1, rho1, m1, impx1, impy1, e1);

		mb = m1;
		impxb = 0.0;
		impyb = impy1;
		eb = e1;
	}

}
void update_boundaries(MPI_Comm grid_comm,
    Field2D<double> &m, Field2D<double> &impx, Field2D<double> &impy, Field2D<double> &e,
	Field2D<int>& klin, int local_Nx, int local_Ny,
    int left, int right, int up, int down,
    int bound_left, int bound_right, int bound_down, int bound_up,
    Params* params) {
	
	double theta = 25.0 * M_PI / 180.0;
	double nx, ny;
	double p_temp1, vx_temp1, vy_temp1, rho_temp1;
	double p_temp2, vx_temp2, vy_temp2, rho_temp2;
    // == ФИЗИЧЕСКИЕ ГРАНИЦЫ ==
    if (left == MPI_PROC_NULL) {
        for (int j = 0; j < local_Ny + 2 * params->fict; ++j) {
            for (int i = 0; i < params->fict; ++i) {
                boundary_cond_x(params->gamma, m(params->fict + i, j), impx(params->fict + i, j), impy(params->fict + i, j), e(params->fict + i, j),
					m(params->fict - i - 1, j), impx(params->fict - i - 1, j), impy(params->fict - i - 1, j), e(params->fict - i - 1, j),
                    bound_left);
            }
        }
    }

    if (right == MPI_PROC_NULL) {
        for (int j = 0; j < local_Ny + 2 * params->fict; ++j) {
            for (int i = 0; i < params->fict; ++i) {
                boundary_cond_x(params->gamma, m(local_Nx + params->fict - i - 1, j), impx(local_Nx + params->fict - i - 1, j), impy(local_Nx + params->fict - i - 1, j),
                    e(local_Nx + params->fict - i - 1, j),
					m(i + local_Nx + params->fict, j), impx(i + local_Nx + params->fict, j), impy(i + local_Nx + params->fict, j), e(i + local_Nx + params->fict, j),
                    bound_right);
            }
        }
    }
    if (down == MPI_PROC_NULL) {
        for (int i = 0; i < local_Nx + 2 * params->fict; ++i) {
            for (int j = 0; j < params->fict; ++j) {
                boundary_cond_y(params->gamma, m(i, params->fict + j), impx(i, params->fict + j), impy(i, params->fict + j), e(i, params->fict + j),
                    m(i, params->fict - j - 1), impx(i, params->fict - j - 1), impy(i, params->fict - j - 1), e(i, params->fict - j - 1),
                    bound_down);
            }
        }
    }
    if (up == MPI_PROC_NULL) {
        for (int i = 0; i < local_Nx + 2 * params->fict; ++i) {
            for (int j = 0; j < params->fict; ++j) {
                boundary_cond_y(params->gamma, m(i, local_Ny - 1 + params->fict - j), impx(i, local_Ny - 1 + params->fict - j), impy(i, local_Ny - 1 + params->fict - j), e(i, local_Ny - 1 + params->fict - j),
                    m(i, local_Ny + params->fict + j), impx(i, local_Ny + params->fict + j), impy(i, local_Ny + params->fict + j), e(i, local_Ny + params->fict + j),
                    bound_up);
            }
        }
    }
	for (int i = params->fict; i < local_Nx + params->fict; ++i) {

		for (int j = params->fict; j < local_Ny + params->fict; ++j) {

			if ((klin(i, j) == 1 && klin(i, j + 1) == 0)) {

				nx = -sin(theta);
				ny = cos(theta);

				for (int k = 0; k < 1; k++) {

					//условие на скорости эквивалентно условию на импульсы
					int i_fluid = i + k;
					int j_fluid = j + k + 1;
					/*std::cout << i << ", " << j - k << std::endl;
					std::cout << i << ", " << j + k + 1 << std::endl;*/
					m(i, j - k) = m(i_fluid, j_fluid);
					e(i, j - k) = e(i_fluid, j_fluid);

					double imp_n = impx(i_fluid, j_fluid) * nx + impy(i_fluid, j_fluid) * ny;
					/*std::cout << v_n << std::endl;*/
					impx(i, j - k) = impx(i_fluid, j_fluid) - 2 * imp_n * nx;
					impy(i, j - k) = impy(i_fluid, j_fluid) - 2 * imp_n * ny;
				}
			}

		}
	}
    // Определение типа данных для столбцов
    MPI_Datatype column_type;
    MPI_Type_vector(local_Ny + 2 * params->fict, 1, local_Nx + 2 * params->fict, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);

    // Уникальные теги для сообщений
    const int TAG_LEFT_TO_RIGHT = 100;
    const int TAG_RIGHT_TO_LEFT = 200;
    const int TAG_DOWN_TO_UP = 300;
    const int TAG_UP_TO_DOWN = 400;

    Field2D<double>* fields[] = { &m, &impx, &impy, &e };
    int k = 1;

    for (auto* f : fields) {
        // Обмен по X (лево-право)
        for (ptrdiff_t i = 0; i < params->fict; ++i) {
            if (right != MPI_PROC_NULL) {
                // Отправляем последние fict расчетных ячеек вправо, принимаем в правые фиктивные ячейки от правого соседа
                MPI_Sendrecv(&(*f)(local_Nx + i, 0), 1, column_type, right, TAG_LEFT_TO_RIGHT * k,
                             &(*f)(local_Nx + params->fict + i, 0), 1, column_type, right, TAG_RIGHT_TO_LEFT * k,
                             grid_comm, MPI_STATUS_IGNORE);
            }
            if (left != MPI_PROC_NULL) {
                // Отправляем первые fict расчетных ячеек влево, принимаем в левые фиктивные ячейки от левого соседа
                MPI_Sendrecv(&(*f)(params->fict + i, 0), 1, column_type, left, TAG_RIGHT_TO_LEFT * k,
                             &(*f)(i, 0), 1, column_type, left, TAG_LEFT_TO_RIGHT * k,
                             grid_comm, MPI_STATUS_IGNORE);
            }
        }

        // Обмен по Y (низ-верх)
        for (ptrdiff_t j = 0; j < params->fict; ++j) {
            if (up != MPI_PROC_NULL) {
                // Отправляем последние fict расчетных строк вверх, принимаем в верхние фиктивные ячейки от верхнего соседа
                MPI_Sendrecv(&(*f)(0, local_Ny + j), local_Nx + 2 * params->fict, MPI_DOUBLE, up, TAG_DOWN_TO_UP * k,
                             &(*f)(0, local_Ny + params->fict + j), local_Nx + 2 * params->fict, MPI_DOUBLE, up, TAG_UP_TO_DOWN * k,
                             grid_comm, MPI_STATUS_IGNORE);
            }
            if (down != MPI_PROC_NULL) {
                // Отправляем первые fict расчетных строк вниз, принимаем в нижние фиктивные ячейки от нижнего соседа
                MPI_Sendrecv(&(*f)(0, params->fict + j), local_Nx + 2 * params->fict, MPI_DOUBLE, down, TAG_UP_TO_DOWN * k,
                             &(*f)(0, j), local_Nx + 2 * params->fict, MPI_DOUBLE, down, TAG_DOWN_TO_UP * k,
                             grid_comm, MPI_STATUS_IGNORE);
            }
        }
        ++k;
    }

    MPI_Type_free(&column_type);
}

void grid(Params* params, vec& x, vec& y, vec& xc, vec& yc, int local_Nx, int local_Ny, double left_b, double down_b)
{
	
	for (int i = 0; i < local_Nx + 2 * params->fict + 1; ++i) {
		x[i] = left_b + (i - params->fict) * params->dx;
	}
	for (int i = 0; i < local_Nx + 2 * params->fict; ++i) {
		xc[i] = 0.5 * (x[i] + x[i + 1]);
	}
	for (int i = 0; i < local_Ny + 2 * params->fict + 1; ++i) {
		y[i] = down_b + (i - params->fict) * params->dy;
	}

	for (int i = 0; i < local_Ny + 2 * params->fict; ++i) {
		yc[i] = 0.5 * (y[i] + y[i + 1]);
	}
}
void init_mach(int local_Nx, int local_Ny, Params* params, vec x, vec y, Field2D<double>& p,
	Field2D<double>& vx, Field2D<double>& vy, Field2D<double>& rho, Field2D<double>& m,
	Field2D<double>& impx, Field2D<double>& impy, Field2D<double>& e, Field2D<int>& klin) {

	double klin_st = 0.5;

	double theta = 25.0 * M_PI / 180.0;

	double mach, p1, vx1, vy1, rho1, p2, vx2, vy2, rho2, y_tg_theta;

	// double mach refraction
	mach = 3.0;
	rho1 = 1.225;

	p1 = 50.0 * 1e3;
	vx1 = mach * std::sqrt(params->gamma * p1 / rho1);
	vy1 = 0.0;

	rho2 = 1.225;
	p2 = 1e3;
	vx2 = 0.0;
	vy2 = 0.0;

	for (int i = 0; i < local_Nx + 2 * params->fict; i++) {
		for (int j = 0; j < local_Ny + 2 * params->fict; j++) {

			p(i, j) = p2;
			vx(i, j) = vx2;
			vy(i, j) = vy2;
			rho(i, j) = rho2;
			klin(i, j) = 0;

			convert_to_conservative(params->gamma, p(i, j), vx(i, j), vy(i, j), rho(i, j), m(i, j), impx(i, j), impy(i, j), e(i, j));
		}
	}

	for (int i = params->fict; i < local_Nx + params->fict; i++) {
		for (int j = params->fict; j < local_Ny + params->fict; j++) {

			double x_otn = x[i] - klin_st;

			y_tg_theta = std::max(0.0, x_otn) * std::tan(theta);

			if ((y[j] - y_tg_theta) < -params->dy && x_otn >= 0) {

				p(i, j) = 2.0 * p1;
				vx(i, j) = 0.0;
				vy(i, j) = 0.0;
				rho(i, j) = 5.0 * rho1;
				klin(i, j) = 1;
			}
			else if ((y[j] - y_tg_theta) >= -params->dy && (y[j] - y_tg_theta) < 0 && x_otn >= 0) {
				klin(i, j) = 1;
				//std::cout << i << ", " << j << std::endl;
			}

			convert_to_conservative(params->gamma, p(i, j), vx(i, j), vy(i, j), rho(i, j), m(i, j), impx(i, j), impy(i, j), e(i, j));

		}
	}

}
void init_krujok(int local_Nx, int local_Ny, int fict, double gamma, vec x, vec y, Field2D<double>& p,
				Field2D<double>& vx, Field2D<double>& vy, Field2D<double>& rho, Field2D<double>& m, 
				Field2D<double>& impx, Field2D<double>& impy, Field2D<double>& e) {
	double R = 0.2;
	double p1, vx1, vy1, rho1, p2, vx2, vy2, rho2, d, gorizont, stena;

	// toro test 1

	rho1 = 1.0;
	vx1 = 0.0;
	p1 = 1.0;
	vy1 = 0.0;
	
	rho2 = 0.125;
	vx2 = 0.0;
	p2 = 0.1;
	vy2 = 0.0;
	
	stena = 0.5;
	
	
	// ne toro test 1

	//p1 = 1.0; // атмосферное давление
	//vx1 = 0.75;
	//vy1 = 0.0;
	//rho1 = 1.0; // плотность воздуха

	//p2 = 1.0; // атмосферное давление
	//vx2 = 0.;
	//vy2 = 0.0;
	//rho2 = 0.125; // плотность воздуха
	
	for (int i = 0; i < local_Nx + 2 * fict; i++) {
		for (int j = 0; j < local_Ny + 2 * fict; j++) {

			d = (x[i] - 0.5) * (x[i] - 0.5) + (y[j] - 0.5) * (y[j] - 0.5);
			
			
			if (x[i] < stena) {
				//d < R*R
			// if (x[i] < stena) {
				
				p(i,j) = p1;
				vx(i, j) = vx1;
				vy(i, j) = vy1;
				rho(i, j) = rho1;
			}
			else {
				p(i, j) = p2;
				vx(i, j) = vx2;
				vy(i, j) = vy2;
				rho(i, j) = rho2;
			}
			convert_to_conservative(gamma, p(i, j), vx(i, j), vy(i, j), rho(i, j), m(i, j), impx(i, j), impy(i, j), e(i, j));
		}
	}
}

double get_dt(Params* params, Field2D<double> p, Field2D<double> rho, Field2D<double> vx, Field2D<double> vy, vec x, vec y, int local_Nx, int local_Ny)
{
	double dt = 10.0e6;
	double c, dt_step;
	for (int i = params->fict; i < local_Nx + params->fict; ++i) {
		for (int j = params->fict; j < local_Ny + params->fict; ++j) {
			c = sqrt(params->gamma * p(i,j) / rho(i, j));
			double denom_x = std::max(std::abs(vx(i, j)) + c, 1e-8);
			double denom_y = std::max(std::abs(vy(i, j)) + c, 1e-8);

			dt_step = std::min(abs(params->CFL * (x[i + 1] - x[i]) / denom_x), abs(params->CFL * (y[j + 1] - y[j]) / denom_y));
			if (dt_step < dt) {
				dt = dt_step;
			}
		}
	}

	return dt;
}

void Godunov_flux_x(double gamma, double p, double vx, double vy, double rho, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
	double m, impx, impy, e;
	//convert_to_conservative(double gamma, double& p, double& vx, double& vy, double& rho, double& m, double& impx, double& impy, double& e)
	convert_to_conservative(gamma, p, vx, vy, rho, m, impx, impy, e);
	Fm = rho * vx;
	Fimpx = Fm * vx + p;
	Fimpy = Fm * vy;
	Fe = (p + e) * vx;
}

void Godunov_flux_y(double gamma, double p, double vx, double vy, double rho, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
	double m, impx, impy, e;
	//convert_to_conservative(double gamma, double& p, double& vx, double& vy, double& rho, double& m, double& impx, double& impy, double& e)
	convert_to_conservative(gamma, p, vx, vy, rho, m, impx, impy, e);
	Fm = rho * vy;
	Fimpx = Fm * vx;
	Fimpy = Fm * vy + p;
	Fe = (p + e) * vy;
}


void WENO_flux(double gamma, double m, double impx, double impy, double e, double c_max, double Fm, double Fimpx, double Fimpy, double Fe,
	double& Fm_minus, double& Fimpx_minus, double& Fimpy_minus, double& Fe_minus, 
	double& Fm_plus, double& Fimpx_plus, double& Fimpy_plus, double& Fe_plus) {

	//расщепление потоков
	Fm_minus = 0.5 * (Fm - c_max * m);
	Fm_plus = 0.5 * (Fm + c_max * m);

	Fimpx_minus = 0.5 * (Fimpx - c_max * impx);
	Fimpx_plus = 0.5 * (Fimpx + c_max * impx);
	
	Fimpy_minus = 0.5 * (Fimpy - c_max * impy);
	Fimpy_plus = 0.5 * (Fimpy + c_max * impy);

	Fe_minus = 0.5 * (Fe - c_max * e);
	Fe_plus = 0.5 * (Fe + c_max * e);

	return;
}
void WENO_reconstructLeft(int fict, int local_Nx, int local_Ny, Field2D<double> v, Field2D<double>& Flow) {

	double eps = 1e-8;
	vec gamma = { 0.3, 0.6, 0.1 };
	vec a(3);
	vec h(3);
	vec beta(3);
	vec omega(3);

	//потоки для разности L


	//задействовано 5 точек
	for (int i = 2; i < local_Nx + 2 * fict - 2; ++i) {
		for (int j = 2; j < local_Ny + 2 * fict - 2; ++j) {
			//шаблоны: для левого

			h[0] = 11.0 / 6.0 * v(i, j) - 7.0 / 6.0 * v(i + 1, j)  + 1.0 / 3.0 * v(i + 2, j) ;

			h[1] = 1.0 / 3.0 * v(i - 1, j)  + 5.0 / 6.0 * v(i, j) - 1.0 / 6.0 * v(i + 1, j) ;

			h[2] = -1.0 / 6.0 * v(i - 2, j)  + 5.0 / 6.0 * v(i - 1, j)  + 1.0 / 3.0 * v(i, j);

			//индикаторы гладкости

			beta[0] = 13.0 / 12.0 * pow(v(i, j)  - 2.0 * v(i + 1, j)  + v(i + 2, j) , 2.0) + 1.0 / 4.0 * pow(3.0 * v(i, j)  - 4.0 * v(i + 1, j)  + v(i + 2, j) , 2.0);

			beta[1] = 13.0 / 12.0 * pow(v(i - 1, j)  - 2.0 * v(i, j)  + v(i + 1, j) , 2.0) + 1.0 / 4.0 * pow(v(i - 1, j)  - v(i + 1, j) , 2.0);

			beta[2] = 13.0 / 12.0 * pow(v(i - 2, j)  - 2.0 * v(i - 1, j)  + v(i, j) , 2.0) + 1.0 / 4.0 * pow(v(i - 2, j)  - 4.0 * v(i - 1, j)  + 3.0 * v(i, j) , 2.0);

			double summ = 0.;

			for (int k = 0; k < 3; ++k) {
				if (beta[k] != 0.) a[k] = gamma[k] / pow(beta[k], 2);

				else  a[k] = gamma[k] / pow(eps, 2.0);

				summ += a[k];
			}

			for (int k = 0; k < 3; ++k) omega[k] = a[k] / summ;

			//считаем левый "поток" 

			Flow(i, j) = std::inner_product(omega.begin(), omega.end(), h.begin(), 0.);
		}
	}
	return;
}
void WENO_reconstructRight(int fict, int local_Nx, int local_Ny, Field2D<double> v, Field2D<double>& Flow) {

	double eps = 1e-8;
	vec gamma = { 0.3, 0.6, 0.1 };
	vec a(3);
	vec h(3);
	vec beta(3);
	vec omega(3);

	//разворачиваем gamma
	reverse(gamma.begin(), gamma.end());

	//потоки для разности L

	//задействовано 5 точек
	for (int i = 2; i < local_Nx + 2 * fict - 2; ++i) {
		for (int j = 2; j < local_Ny + 2 * fict - 2; ++j) {
			//шаблоны: для правого

			h[0] = 1.0 / 3.0 * v(i, j)  + 5.0 / 6.0 * v(i + 1, j) - 1.0 / 6.0 * v(i + 2, j);

			h[1] = -1.0 / 6.0 * v(i - 1, j) + 5.0 / 6.0 * v(i, j)  + 1.0 / 3.0 * v(i + 1, j);

			h[2] = 1.0 / 3.0 * v(i - 2, j) - 7.0 / 6.0 * v(i - 1, j) + 11.0 / 6.0 * v(i, j) ;

			//индикаторы гладкости
			beta[0] = 13.0 / 12.0 * pow(v(i, j)  - 2.0 * v(i + 1, j) + v(i + 2, j), 2.0) + 1.0 / 4.0 * pow(3.0 * v(i, j)  - 4.0 * v(i + 1, j) + v(i + 2, j), 2.0);

			beta[1] = 13.0 / 12.0 * pow(v(i - 1, j) - 2.0 * v(i, j)  + v(i + 1, j), 2.0) + 1.0 / 4.0 * pow(v(i - 1, j) - v(i + 1, j), 2.0);

			beta[2] = 13.0 / 12.0 * pow(v(i - 2, j) - 2.0 * v(i - 1, j) + v(i, j) , 2.0) + 1.0 / 4.0 * pow(v(i - 2, j) - 4.0 * v(i - 1, j) + 3.0 * v(i, j) , 2.0);

			double summ = 0.;

			for (int k = 0; k < 3; ++k) {
				if (beta[k] != 0.) a[k] = gamma[k] / pow(beta[k], 2);

				else  a[k] = gamma[k] / pow(eps, 2.0);

				summ += a[k];
			}

			for (int k = 0; k < 3; ++k) omega[k] = a[k] / summ;

			//считаем правый "поток" 

			Flow(i, j) = inner_product(omega.begin(), omega.end(), h.begin(), 0.);
		}
	}
	return;
}
void WENO_reconstructDown(int fict, int local_Nx, int local_Ny, Field2D<double> v, Field2D<double>& Flow) {

	double eps = 1e-8;
	vec gamma = { 0.3, 0.6, 0.1 };
	vec a(3);
	vec h(3);
	vec beta(3);
	vec omega(3);

	//потоки для разности L


	//задействовано 5 точек
	for (int i = 2; i < local_Nx + 2 * fict - 2; ++i) {
		for (int j = 2; j < local_Ny + 2 * fict - 2; ++j) {
			//шаблоны: для левого

			h[0] = 11.0 / 6.0 * v(i, j) - 7.0 / 6.0 * v(i, j + 1) + 1.0 / 3.0 * v(i, j + 2);

			h[1] = 1.0 / 3.0 * v(i, j - 1) + 5.0 / 6.0 * v(i, j) - 1.0 / 6.0 * v(i, j + 1);

			h[2] = -1.0 / 6.0 * v(i, j - 2) + 5.0 / 6.0 * v(i, j - 1) + 1.0 / 3.0 * v(i, j);

			//индикаторы гладкости

			beta[0] = 13.0 / 12.0 * pow(v(i, j) - 2.0 * v(i, j + 1) + v(i, j + 2), 2.0) + 1.0 / 4.0 * pow(3.0 * v(i, j) - 4.0 * v(i, j + 1) + v(i, j + 2), 2.0);

			beta[1] = 13.0 / 12.0 * pow(v(i, j - 1) - 2.0 * v(i, j) + v(i, j + 1), 2.0) + 1.0 / 4.0 * pow(v(i, j - 1) - v(i, j + 1), 2.0);

			beta[2] = 13.0 / 12.0 * pow(v(i, j - 2) - 2.0 * v(i, j - 1) + v(i, j), 2.0) + 1.0 / 4.0 * pow(v(i, j - 2) - 4.0 * v(i, j - 1) + 3.0 * v(i, j), 2.0);

			double summ = 0.;

			for (int k = 0; k < 3; ++k) {
				if (beta[k] != 0.) a[k] = gamma[k] / pow(beta[k], 2);

				else  a[k] = gamma[k] / pow(eps, 2.0);

				summ += a[k];
			}

			for (int k = 0; k < 3; ++k) omega[k] = a[k] / summ;

			//считаем левый "поток" 

			Flow(i, j) = std::inner_product(omega.begin(), omega.end(), h.begin(), 0.);
		}
	}
	return;
}
void WENO_reconstructUp(int fict, int local_Nx, int local_Ny, Field2D<double> v, Field2D<double>& Flow) {

	double eps = 1e-8;
	vec gamma = { 0.3, 0.6, 0.1 };
	vec a(3);
	vec h(3);
	vec beta(3);
	vec omega(3);

	//разворачиваем gamma
	reverse(gamma.begin(), gamma.end());

	//потоки для разности L

	//задействовано 5 точек
	for (int i = 2; i < local_Nx + 2 * fict - 2; ++i) {
		for (int j = 2; j < local_Ny + 2 * fict - 2; ++j) {
			//шаблоны: для правого

			h[0] = 1.0 / 3.0 * v(i, j) + 5.0 / 6.0 * v(i, j + 1) - 1.0 / 6.0 * v(i, j + 2);

			h[1] = -1.0 / 6.0 * v(i, j - 1) + 5.0 / 6.0 * v(i, j) + 1.0 / 3.0 * v(i, j + 1);

			h[2] = 1.0 / 3.0 * v(i, j - 2) - 7.0 / 6.0 * v(i, j - 1) + 11.0 / 6.0 * v(i, j);

			//индикаторы гладкости
			beta[0] = 13.0 / 12.0 * pow(v(i, j) - 2.0 * v(i, j + 1) + v(i, j + 2), 2.0) + 1.0 / 4.0 * pow(3.0 * v(i, j) - 4.0 * v(i, j + 1) + v(i, j + 2), 2.0);

			beta[1] = 13.0 / 12.0 * pow(v(i, j - 1) - 2.0 * v(i, j) + v(i, j + 1), 2.0) + 1.0 / 4.0 * pow(v(i, j - 1) - v(i, j + 1), 2.0);

			beta[2] = 13.0 / 12.0 * pow(v(i, j - 2) - 2.0 * v(i, j - 1) + v(i, j), 2.0) + 1.0 / 4.0 * pow(v(i, j - 2) - 4.0 * v(i, j - 1) + 3.0 * v(i, j), 2.0);

			double summ = 0.;

			for (int k = 0; k < 3; ++k) {
				if (beta[k] != 0.) a[k] = gamma[k] / pow(beta[k], 2);

				else  a[k] = gamma[k] / pow(eps, 2.0);

				summ += a[k];
			}

			for (int k = 0; k < 3; ++k) omega[k] = a[k] / summ;

			//считаем правый "поток" 

			Flow(i, j) = inner_product(omega.begin(), omega.end(), h.begin(), 0.);
		}
	}
	return;
}
double Sound_velocity(double gamma, double p_sound_velocity, double rho_sound_velocity) {
	return std::sqrt(gamma * p_sound_velocity / rho_sound_velocity);
}
void Artificial_viscosity(int viscosity_type, Field2D<double>& p, Field2D<double>& vx, Field2D<double>& vy,
	Field2D<double>& rho, Field2D<double>& m_next, Field2D<double>& impx_next, Field2D<double>& impy_next, Field2D<double>& e_next,
	int local_Nx, int local_Ny, Params* params) {

	double mu_1 = 0.1;
	double mu_2 = 0.2;
	double mu_3 = 0.2;

	for (int i = params->fict; i < local_Nx + params->fict; ++i) {
		for (int j = params->fict; j < local_Ny + params->fict; ++j) {
			double wx = 0.0, wy = 0.0, p_next = 0.0;
			//физическая, на участках уменьшения скорости
			if (viscosity_type == 1) {
				if (vx(i + 1, j) - vx(i, j) < 0) {
					wx = -mu_1 * rho(i, j) * (vx(i + 1, j) - vx(i, j));
				}
				if (vy(i, j + 1) - vy(i, j) < 0) {
					wy = -mu_1 * rho(i, j) * (vy(i, j + 1) - vy(i, j));
				}
				p(i, j) += wx + wy;
			}
			//линейная - первая степень по x, по vx, но со скоростью звука
			else if (viscosity_type == 2) {
				double c = sqrt(params->gamma * p(i, j) / rho(i, j));

				if (vx(i + 1, j) - vx(i, j) < 0) {
					wx = -mu_2 * rho(i, j) * c * (params->dx) * (vx(i + 1, j) - vx(i, j));
				}
				if (vy(i, j + 1) - vy(i, j) < 0) {
					wy = -mu_2 * rho(i, j) * c * (params->dy) * (vy(i, j + 1) - vy(i, j));
				}
				p(i, j) += wx + wy;
			}
			// Неймана-Рихтмайера - квадрат скорости и шага по простр
			else if (viscosity_type == 3) {
				if (vx(i + 1, j) - vx(i, j) < 0) {
					wx = -mu_3 * rho(i, j) * abs(params->dx) * (params->dx) * abs((vx(i + 1, j) - vx(i, j))) * (vx(i + 1, j) - vx(i, j));
				}
				if (vy(i, j + 1) - vy(i, j) < 0) {
					wy = -mu_3 * rho(i, j) * abs(params->dy) * (params->dy) * abs((vy(i, j + 1) - vy(i, j))) * (vy(i, j + 1) - vy(i, j));
				}
				p(i, j) += wx + wy;
			}
			// Куропатенко - нет констант
			else if (viscosity_type == 4) {
				double c = sqrt(params->gamma * p(i, j) / rho(i, j));
				if (vx(i + 1, j) - vx(i, j) < 0) {

					wx = (params->gamma + 1) / 4 * rho(i, j) * pow(params->dx * (vx(i + 1, j) - vx(i, j)), 2);

					wx = wx + sqrt(pow(wx, 2) + (params->gamma + 1) / 4 * pow(rho(i, j) * c * params->dx * (vx(i + 1, j) - vx(i, j)), 2));
				}
				if (vy(i, j + 1) - vy(i, j) < 0) {
					wy = (params->gamma + 1) / 4 * rho(i, j) * pow(params->dy * (vy(i, j + 1) - vy(i, j)), 2);

					wy = wy + sqrt(pow(wy, 2) + (params->gamma + 1) / 4 * pow(rho(i, j) * c * params->dy * (vy(i, j + 1) - vy(i, j)), 2));
				}
			}

			e_next(i, j) += (wx + wy) / (params->gamma - 1.0);
		}
	}

}
void GK_2d(int argc, char** argv)
{
	int rank, size;
	MPI_Init(&argc, &argv);

	MPI_Status Status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	std::vector<int> dims(2, 0);
	//MPI_Dims_create(
	//	size, // Общее число процессов
	//	2, // Количество измерений
	//	dims.data() // Результат разбиения
	//);

	dims = delimetrs(size);
	std::cout << dims[0] << ", " << dims[1] << std::endl;
	int periods[2] = { 0, 0 };  // без периодических границ у Полины было так
	int reorder = 1;          // позволить MPI изменить порядок номеров процессов для эффективности
	
	MPI_Comm grid_comm;
	MPI_Cart_create(
		MPI_COMM_WORLD,
		2,        // размерность сетки (x, y)
		dims.data(),     // разбиение от Dims_create
		periods,  // граничные условия
		reorder,  // можно ли менять ранги процессов
		&grid_comm
	);

	//проверка
	if (grid_comm != MPI_COMM_NULL) {
		std::vector<int> coords(2, 0);
		MPI_Cart_coords(grid_comm, rank, 2, coords.data());
		// Дальнейшая логика только для процессов в grid_comm
	}
	
	// Узнаем координаты текущего процесса в сетке
	std::vector<int> coords(2, 0);
	MPI_Cart_coords(grid_comm, rank, 2, coords.data());

	int left, right, up, down;
	MPI_Cart_shift(grid_comm, 0, 1, &left, &right); // 0 = направление X
	MPI_Cart_shift(grid_comm, 1, 1, &down, & up);    // 1 = направление Y
	//std::cout << "Rank " << rank << ", left = " << left << ", MPI_PROC_NULL = " << MPI_PROC_NULL << ", (equal?) " << (left == MPI_PROC_NULL) << std::endl;
	//std::cout << "Rank " << rank << ", right = " << right << ", MPI_PROC_NULL = " << MPI_PROC_NULL << ", (equal?) " << (right == MPI_PROC_NULL) << std::endl;
	
	//тип указал 2 - для ВЕНО. Подправил парамс, чтобы при type >= 2 fict = 3
	Params* params = new Params(std::stol(argv[1]), std::stol(argv[2]), 2, 0.15, 4.0, 2.0, 1.0, 1.4, 2.0);
	std::cout<< params->Nx << ", " << params->Ny << ", " << params->T <<std::endl;

	// Базовый размер блока
	int base_x = params->Nx / dims[0];
	int base_y = params->Ny / dims[1];
	// Остатки
	int extra_x = params->Nx % dims[0];
	int extra_y = params->Ny % dims[1];
	// Стартовые индексы
	int start_x = coords[0] * base_x + (coords[0] < extra_x ? coords[0] : extra_x);
	int start_y = coords[1] * base_y + (coords[1] < extra_y ? coords[1] : extra_y);
	// Конечные индексы
	int end_x = start_x + base_x + (coords[0] < extra_x ? 1 : 0);
	int end_y = start_y + base_y + (coords[1] < extra_y ? 1 : 0);

	int local_Nx = end_x - start_x;
	int local_Ny = end_y - start_y;

	// пересмотреть где какие циферки (+- fict и сколько) задавать
	
	int xsize = local_Nx + 2 * params->fict, ysize = local_Ny + 2 * params->fict;
	
	vec xc(xsize);
	vec x(xsize + 1);
	vec yc(ysize);
	vec y(ysize + 1);

	//+ 2 для оптимизации обмена между процессами
	

	Field2D<double> p(xsize, ysize), vx(xsize, ysize), vy(xsize, ysize), rho(xsize, ysize);
	Field2D<double> m(xsize, ysize), impx(xsize, ysize), impy(xsize, ysize), e(xsize, ysize);

	Field2D<double> Fm_plus_x(xsize, ysize), Fimpx_plus_x(xsize, ysize),
					Fimpy_plus_x(xsize, ysize), Fe_plus_x(xsize, ysize);
	Field2D<double>Fm_minus_x(xsize, ysize), Fimpx_minus_x(xsize, ysize),
					Fimpy_minus_x(xsize, ysize), Fe_minus_x(xsize, ysize);

	Field2D<double> Fm_plus_y(xsize, ysize), Fimpx_plus_y(xsize, ysize),
		Fimpy_plus_y(xsize, ysize), Fe_plus_y(xsize, ysize);
	Field2D<double>Fm_minus_y(xsize, ysize), Fimpx_minus_y(xsize, ysize),
		Fimpy_minus_y(xsize, ysize), Fe_minus_y(xsize, ysize);

	Field2D<double> Flowm_l(xsize, ysize), Flowm_r(xsize, ysize),
		Flowm_d(xsize, ysize), Flowm_u(xsize, ysize);
	Field2D<double>Flowimpx_l(xsize, ysize), Flowimpx_r(xsize, ysize),
		Flowimpx_d(xsize, ysize), Flowimpx_u(xsize, ysize);
	Field2D<double>Flowimpy_l(xsize, ysize), Flowimpy_r(xsize, ysize),
		Flowimpy_d(xsize, ysize), Flowimpy_u(xsize, ysize);
	Field2D<double>Flowe_l(xsize, ysize), Flowe_r(xsize, ysize),
		Flowe_d(xsize, ysize), Flowe_u(xsize, ysize);

	Field2D<double>Lm(xsize, ysize), Limpx(xsize, ysize),
		Limpy(xsize, ysize), Le(xsize, ysize);
	double Fm, Fimpx, Fimpy, Fe;

	// задание границ клина
	Field2D<int> klin(local_Nx + 2 * params->fict, local_Ny + 2 * params->fict);
	std::vector<vec> rk = { {1.0,    0.0,       0.0,       0.0},
									  {1.0,    0.0,       0.0,       1.0},
									  {0.75, 0.25,       0.0,    0.25},
									  {1.0 / 3.0,    0.0, 2.0 / 3.0, 2.0 / 3.0} };

	int step = 0, max_step = 10000;// 201;
	double time = 0.0;
	double dt_temp;

	int out = 100;
	int bound_left = 2, bound_right = 1, bound_down = 0, bound_up = 0;
	int viscosity_type = 0;

	int left_b = 0, down_b = 0;
	//MPI_Exscan(&local_Nx, &left_b, 1, MPI_INT, MPI_SUM, grid_comm);
	//MPI_Exscan(&local_Ny, &down_b, 1, MPI_INT, MPI_SUM, grid_comm);

	// Создаем подкоммуникаторы для строк и столбцов
    	MPI_Comm row_comm, col_comm;

    	// 1. Коммуникатор для строк (фиксируем Y, меняется X)
    	MPI_Comm_split(grid_comm, coords[1], coords[0], &row_comm);

    	// 2. Коммуникатор для столбцов (фиксируем X, меняется Y)
    	MPI_Comm_split(grid_comm, coords[0], coords[1], &col_comm);

   	 // Вычисляем left_b (только внутри строки)
 	MPI_Exscan(&local_Nx, &left_b, 1, MPI_INT, MPI_SUM, row_comm);
  	if (coords[0] == 0) left_b = 0; // Первый столбец начинается с 0

    	// Вычисляем down_b (только внутри столбца)
    	MPI_Exscan(&local_Ny, &down_b, 1, MPI_INT, MPI_SUM, col_comm);
    	if (coords[1] == 0) down_b = 0; // Первая строка начинается с 0

    	// Освобождаем коммуникаторы
    	MPI_Comm_free(&row_comm);
    	MPI_Comm_free(&col_comm);

		//std::cout << "splits ok" << std::endl;
	grid(params, x, y, xc, yc, local_Nx, local_Ny, left_b * (params->dx), down_b * (params->dy));
	//std::cout << "grid ok" << std::endl;
	init_mach(local_Nx, local_Ny, params, x, y, p, vx, vy, rho, m, impx, impy, e, klin);
		
	//std::cout << "init ok" << std::endl;
	
	//std::cout << rank << " " << rho(0, 5) << " " << rho(xsize - 1, 5) << std::endl;
	//std::cout << rank << " " << rho(1, 5) << " " << rho(xsize - 2, 5) << std::endl;
	data_to_file(local_Nx, local_Ny, params->fict, params->gamma, 0.0, x, y, p, vx, vy, rho, rank, dims, (argv[3]));
	
	//exit(1);
	//std::cout << "data ok" << std::endl;
	//временные массивы для РК
	std::vector<Field2D<double>> m_temp(4, Field2D<double>(xsize, ysize));
	std::vector<Field2D<double>> impx_temp(4, Field2D<double>(xsize, ysize));
	std::vector<Field2D<double>> impy_temp(4, Field2D<double>(xsize, ysize));
	std::vector<Field2D<double>> e_temp(4, Field2D<double>(xsize, ysize));
	
	while (time < params->T && step < max_step)
	{	
		/*if (step % 100 == 0) {
			log_ghost_cells(rank, step, rho, local_Nx, local_Ny, params->fict);
		}*/
		update_boundaries(grid_comm, m, impx, impy, e, klin, local_Nx, local_Ny,
			left, right, up, down,
			bound_left, bound_right, bound_down, bound_up,
			params);
		
		m_temp[0] = m;
		impx_temp[0] = impx;
		impy_temp[0] = impy;
		e_temp[0] = e;

		//std::cout << "update ok" << std::endl;
		dt_temp = get_dt(params, p, rho, vx, vy, x, y, local_Nx, local_Ny);
		
		double dt = 0.0;

		MPI_Allreduce(&dt_temp, &dt, 1, MPI_DOUBLE, MPI_MIN, grid_comm);

		time += dt;
		//std::cout << dt << std::endl;
		/*if (step == 72) {
			std::cout << step << std::endl;
		}*/
		double A_temp = 0., A = 0.;
		for (int i = params->fict; i < local_Nx + params->fict; ++i) {
			for (int j = params->fict; j < local_Ny + params->fict; ++j) {
				// максимальная хар. скорость
				double c = Sound_velocity(params->gamma, (params->gamma - 1.0) * (e(i, j) - 0.5 * rho(i, j) * (vx(i, j) * vx(i, j) +
					vy(i, j) * vy(i, j))), rho(i, j));
				A_temp = std::max(A_temp, c);
			}
		}

		MPI_Allreduce(&A_temp, &A, 1, MPI_DOUBLE, MPI_MAX, grid_comm);

		/*   START CALC   */
		for (int s = 1; s < 4; ++s) {

			for (int i = 0; i < local_Nx + 2 * params->fict; i++) {

				for (int j = 0; j < local_Ny + 2 * params->fict; j++) {

					double p_temp, vx_temp, vy_temp, rho_temp;
					convert_from_conservative(params->gamma, p_temp, vx_temp, vy_temp, rho_temp,
						m_temp[s - 1](i,j), impx_temp[s - 1](i, j), impy_temp[s - 1](i, j), e_temp[s - 1](i, j));


					//расщепление по х
					Godunov_flux_x(params->gamma, p_temp, vx_temp, vy_temp, rho_temp, Fm, Fimpx, Fimpy, Fe);

					WENO_flux(params->gamma, m_temp[s - 1](i, j), impx_temp[s - 1](i, j), impy_temp[s - 1](i, j), e_temp[s - 1](i, j),
							A, Fm, Fimpx, Fimpy, Fe, Fm_minus_x(i, j), Fimpx_minus_x(i, j), Fimpy_minus_x(i, j), Fe_minus_x(i, j),
							Fm_plus_x(i, j), Fimpx_plus_x(i, j), Fimpy_plus_x(i, j), Fe_plus_x(i, j));

					//расщепление по у
					Godunov_flux_y(params->gamma, p_temp, vx_temp, vy_temp, rho_temp, Fm, Fimpx, Fimpy, Fe);

					WENO_flux(params->gamma, m_temp[s - 1](i, j), impx_temp[s - 1](i, j), impy_temp[s - 1](i, j), e_temp[s - 1](i, j),
						A, Fm, Fimpx, Fimpy, Fe, Fm_minus_y(i, j), Fimpx_minus_y(i, j), Fimpy_minus_y(i, j), Fe_minus_y(i, j),
						Fm_plus_y(i, j), Fimpx_plus_y(i, j), Fimpy_plus_y(i, j), Fe_plus_y(i, j));

				}
			}
			WENO_reconstructLeft(params->fict, local_Nx, local_Ny, Fm_minus_x, Flowm_l);
			WENO_reconstructRight(params->fict, local_Nx, local_Ny, Fm_plus_x, Flowm_r);
			WENO_reconstructDown(params->fict, local_Nx, local_Ny, Fm_minus_y, Flowm_d);
			WENO_reconstructUp(params->fict, local_Nx, local_Ny, Fm_plus_y, Flowm_u);

			WENO_reconstructLeft(params->fict, local_Nx, local_Ny, Fimpx_minus_x, Flowimpx_l);
			WENO_reconstructRight(params->fict, local_Nx, local_Ny, Fimpx_plus_x, Flowimpx_r);
			WENO_reconstructDown(params->fict, local_Nx, local_Ny, Fimpx_minus_y, Flowimpx_d);
			WENO_reconstructUp(params->fict, local_Nx, local_Ny, Fimpx_plus_y, Flowimpx_u);

			WENO_reconstructLeft(params->fict, local_Nx, local_Ny, Fimpy_minus_x, Flowimpy_l);
			WENO_reconstructRight(params->fict, local_Nx, local_Ny, Fimpy_plus_x, Flowimpy_r);
			WENO_reconstructDown(params->fict, local_Nx, local_Ny, Fimpy_minus_y, Flowimpy_d);
			WENO_reconstructUp(params->fict, local_Nx, local_Ny, Fimpy_plus_y, Flowimpy_u);

			WENO_reconstructLeft(params->fict, local_Nx, local_Ny, Fe_minus_x, Flowe_l);
			WENO_reconstructRight(params->fict, local_Nx, local_Ny, Fe_plus_x, Flowe_r);
			WENO_reconstructDown(params->fict, local_Nx, local_Ny, Fe_minus_y, Flowe_d);
			WENO_reconstructUp(params->fict, local_Nx, local_Ny, Fe_plus_y, Flowe_u);

			//считаем смещение для РК по x
			for (int i = 1; i < local_Nx + 2 * params->fict - 1; ++i) {
				for (int j = 1; j < local_Ny + 2 * params->fict - 1; ++j) {

					double dFmdx = ((Flowm_r(i, j) + Flowm_l(i + 1, j)) - (Flowm_r(i - 1, j) + Flowm_l(i, j))) / params->dx;
					double dGmdy = ((Flowm_u(i, j) + Flowm_d(i, j + 1)) - (Flowm_u(i, j - 1) + Flowm_d(i, j))) / params->dy;
					
					double dFimpxdx = ((Flowimpx_r(i, j) + Flowimpx_l(i + 1, j)) - (Flowimpx_r(i - 1, j) + Flowimpx_l(i, j))) / params->dx;
					double dGimpxdy = ((Flowimpx_u(i, j) + Flowimpx_d(i, j + 1)) - (Flowimpx_u(i, j - 1) + Flowimpx_d(i, j))) / params->dy;
					
					double dFimpydx = ((Flowimpy_r(i, j) + Flowimpy_l(i + 1, j)) - (Flowimpy_r(i - 1, j) + Flowimpy_l(i, j))) / params->dx;
					double dGimpydy = ((Flowimpy_u(i, j) + Flowimpy_d(i, j + 1)) - (Flowimpy_u(i, j - 1) + Flowimpy_d(i, j))) / params->dy;
					
					double dFedx = ((Flowe_r(i, j) + Flowe_l(i + 1, j)) - (Flowe_r(i - 1, j) + Flowe_l(i, j))) / params->dx;
					double dGedy = ((Flowe_u(i, j) + Flowe_d(i, j + 1)) - (Flowe_u(i, j - 1) + Flowe_d(i, j))) / params->dy;
					
					Lm(i, j) = dFmdx + dGmdy;

					Limpx(i, j) = dFimpxdx + dGimpxdy;

					Limpy(i, j) = dFimpydx + dGimpydy;

					Le(i, j) = dFedx + dGedy;

					m_temp[s](i, j) = rk[s][0] * m_temp[0](i, j)
						+ rk[s][1] * m_temp[1](i, j)
						+ rk[s][2] * m_temp[2](i, j)
						- rk[s][3] * dt * Lm(i, j);

					impx_temp[s](i, j) = rk[s][0] * impx_temp[0](i, j)
						+ rk[s][1] * impx_temp[1](i, j)
						+ rk[s][2] * impx_temp[2](i, j)
						- rk[s][3] * dt * Limpx(i, j);

					impy_temp[s](i, j) = rk[s][0] * impy_temp[0](i, j)
						+ rk[s][1] * impy_temp[1](i, j)
						+ rk[s][2] * impy_temp[2](i, j)
						- rk[s][3] * dt * Limpy(i, j);

					e_temp[s](i, j) = rk[s][0] * e_temp[0](i, j)
						+ rk[s][1] * e_temp[1](i, j)
						+ rk[s][2] * e_temp[2](i, j)
						- rk[s][3] * dt * Le(i, j);

					

					if (std::isnan(m_temp[s](i, j))) {
						std::cerr << "‼️ NaN in m_temp[" << s << "] at (" << i << "," << j << "), step = " << step << std::endl;
					}
					if (std::isnan(impx_temp[s](i, j))) {
						std::cerr << "‼️ NaN in impx_temp[" << s << "] at (" << i << "," << j << "), step = " << step << std::endl;
					}
					if (std::isnan(impy_temp[s](i, j))) {
						std::cerr << "‼️ NaN in impy_temp[" << s << "] at (" << i << "," << j << "), step = " << step << std::endl;
					}
					if (std::isnan(e_temp[s](i, j))) {
						std::cerr << "‼️ NaN in e_temp[" << s << "] at (" << i << "," << j << "), step = " << step << std::endl;
					}
				}
			}

			update_boundaries(grid_comm, m_temp[s], impx_temp[s], impy_temp[s], e_temp[s],
				klin, local_Nx, local_Ny, left, right, up, down, bound_left, bound_right, bound_down, bound_up, params);
		}
		
		/*    END CALC    */
		if (viscosity_type != 0)
		{
			Artificial_viscosity(viscosity_type, p, vx, vy, rho, m_temp[3], impx_temp[3], impy_temp[3], e_temp[3], local_Nx, local_Ny, params);
		}

		// update parameters
		for (int i = params->fict; i < local_Nx + params->fict; ++i)
		{
			for (int j = params->fict; j < local_Ny + params->fict; ++j)
			{
				if (klin(i, j) == 0) {
					
					m(i, j) = m_temp[3](i, j);
					impx(i, j) = impx_temp[3](i, j);
					impy(i, j) = impy_temp[3](i, j);
					e(i, j) = e_temp[3](i, j);

					convert_from_conservative(params->gamma, p(i, j), vx(i, j), vy(i, j), rho(i, j), m(i, j), impx(i, j), impy(i, j), e(i, j));
				}
			}
		}
		

		if (step % out == 0)
			data_to_file(local_Nx, local_Ny, params->fict, params->gamma, time, x, y, p, vx, vy, rho, rank, dims, (argv[3]));
		step += 1;
		if (rank == 0 && !(step % out)) { 	std::cout << step << std::endl;}

		

	}
	MPI_Finalize();
	return;
}


int main(int argc, char** argv)
{
	GK_2d(argc, argv);

	return 0;
}