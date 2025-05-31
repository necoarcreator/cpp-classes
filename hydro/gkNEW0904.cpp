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
#include <fstream>
#include <iostream>
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
		else this->fict = 2;
		
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
	p = (gamma - 1.0) * (e - 0.5 * (pow(impx, 2.0) + pow(impy, 2.0)) / m);
	vx = impx / m;
	vy = impy / m;
	rho = m;
}
// to conservative
void convert_to_conservative(double gamma, double& p, double& vx, double& vy, double& rho, double& m, double& impx, double& impy, double& e)
{
	m = rho;
	impx = rho * vx;
	impy = rho * vy;
	e = 0.5 * rho * (pow(vx, 2.0) + pow(vy, 2.0)) + p / (gamma - 1.0);
}


void boundary_cond_x(double gamma, double p, double vx, double vy, double rho, double& pb, double& vxb, double& vyb, double& rhob, int mode)
{
	double mach, rho1, p1, vx1, vy1;

	mach = 3.0;
	rho1 = 1.225;

	p1 = 50.0 * 1e3;
	vx1 = mach * std::sqrt(gamma * p1 / rho1);
	vy1 = 0.0;

	// wall
	if (mode == 0) {
		pb = p;
		vxb = -vx;
		vyb = vy;
		rhob = rho;
	}

	// free flux
	else if (mode == 1) {
		pb = p;
		vxb = vx;
		vyb = vy;
		rhob = rho;
	}

	//mach params
	else {
		pb = p1;
		rhob = rho1;
		vxb = vx1;
		vyb = vy1;
	}

}

void boundary_cond_y(double gamma, double p, double vx, double vy, double rho, double& pb, double& vxb, double& vyb, double& rhob, int mode)
{
	// wall
	if (mode == 0) {
		pb = p;
		vxb = vx;
		vyb = -vy;
		rhob = rho;
	}

	// free flux
	else if (mode == 1) {
		pb = p;
		vxb = vx;
		vyb = vy;
		rhob = rho;
	}

	//mach params
	else {
		pb = 50 * 1e3;
		rhob = 1.225;
		vxb = 0.0;
		vyb = 3.0 * std::sqrt(gamma * pb / rhob);
	}
}

void update_boundaries(MPI_Comm grid_comm, vec x, vec y, 
    Field2D<double> &p, Field2D<double> &vx, Field2D<double> &vy, Field2D<double> &rho,
    Field2D<double> &m, Field2D<double> &impx, Field2D<double> &impy, Field2D<double> &e,
	Field2D<int>& klin, int local_Nx, int local_Ny,
    int left, int right, int up, int down,
    int bound_left, int bound_right, int bound_down, int bound_up,
    Params* params) {

	double theta = 25.0 * M_PI / 180.0;
	double nx, ny;
			
    // == ФИЗИЧЕСКИЕ ГРАНИЦЫ ==
    if (left == MPI_PROC_NULL) {
        for (int j = 0; j < local_Ny + 2 * params->fict; ++j) {
    		for (int i = 0; i < params->fict; ++i) {
        		boundary_cond_x(params->gamma, p(params->fict + i, j), vx(params->fict + i, j), vy(params->fict + i, j), rho(params->fict + i, j),
            		p(params->fict - i - 1, j), vx(params->fict - i - 1, j), vy(params->fict - i - 1, j), rho(params->fict - i - 1, j),
            		bound_left);
    		}
    	}
    }
    if (right == MPI_PROC_NULL) {
        for (int j = 0; j < local_Ny + 2 * params->fict; ++j) {
            for (int i = 0; i < params->fict; ++i) {
                boundary_cond_x(params->gamma, p(local_Nx + params->fict - i - 1, j), vx(local_Nx + params->fict - i - 1, j), vy(local_Nx + params->fict - i - 1, j),
                    rho(local_Nx + params->fict - i - 1, j),
                    p(i + local_Nx + params->fict, j), vx(i + local_Nx + params->fict, j), vy(i + local_Nx + params->fict, j), rho(i + local_Nx + params->fict, j),
                    bound_right);
            }
        }
    }
    if (down == MPI_PROC_NULL) {
        for (int i = 0; i < local_Nx + 2 * params->fict; ++i) {
            for (int j = 0; j < params->fict; ++j) {
                boundary_cond_y(params->gamma, p(i, params->fict + j), vx(i, params->fict + j), vy(i, params->fict + j), rho(i, params->fict + j),
                    p(i, params->fict - j - 1), vx(i, params->fict - j - 1), vy(i, params->fict - j - 1), rho(i, params->fict - j - 1),
                    bound_down);
            }
        }
    }
    if (up == MPI_PROC_NULL) {
        for (int i = 0; i < local_Nx + 2 * params->fict; ++i) {
            for (int j = 0; j < params->fict; ++j) {
                boundary_cond_y(params->gamma, p(i, local_Ny - 1 + params->fict - j), vx(i, local_Ny - 1 + params->fict - j), vy(i, local_Ny - 1 + params->fict - j), rho(i, local_Ny - 1 + params->fict - j),
                    p(i, local_Ny + params->fict + j), vx(i, local_Ny + params->fict + j), vy(i, local_Ny + params->fict + j), rho(i, local_Ny + params->fict + j),
                    bound_up);
            }
        }
    }
    for (int i = params->fict; i < local_Nx + params->fict; ++i) {
    	for (int j = params->fict; j < local_Ny + params->fict; ++j) {
			
			if ((klin(i, j) == 1 && klin(i, j + 1) == 0) || (klin(i, j) == 1 && klin(i - 1, j) == 0)) {

				nx = -sin(theta);
				ny = cos(theta); 

				for (int k = 0; k < 1; k++) {

					int i_fluid = i + k;
					int j_fluid = j + k + 1;
					/*std::cout << i << ", " << j - k << std::endl;
					std::cout << i << ", " << j + k + 1 << std::endl;*/
					p(i, j - k) = p(i_fluid, j_fluid);
					rho(i, j - k) = rho(i_fluid, j_fluid);

					double v_n = vx(i_fluid, j_fluid) * nx + vy(i_fluid, j_fluid) * ny;
					/*std::cout << v_n << std::endl;*/
					vx(i, j - k) = vx(i_fluid, j_fluid) - 2 * v_n * nx;
					vy(i, j - k) = vy(i_fluid, j_fluid) - 2 * v_n * ny;
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

    Field2D<double>* fields[] = { &p, &vx, &vy, &rho };
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

    // Обновление консервативных переменных только в фиктивных ячейках
    for (int i = 0; i < local_Nx + 2 * params->fict; ++i) {
        for (int j = 0; j < local_Ny + 2 * params->fict; ++j) {
            //if (i < params->fict || i >= local_Nx + params->fict || j < params->fict || j >= local_Ny + params->fict || (klin(i, j) == 0 && klin(i, j - 1) == 1)) {
                convert_to_conservative(params->gamma, p(i, j), vx(i, j), vy(i, j), rho(i, j),
                    m(i, j), impx(i, j), impy(i, j), e(i, j));
           // }
        }
    }


}
//*/
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
	p2 = 50.0 * 1e3;
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
			else if ((y[j] - y_tg_theta) >= -params->dy && (y[j] - y_tg_theta)  < 0 && x_otn >= 0) {
				klin(i, j) = 1;
				//std::cout << i << ", " << j << std::endl;
			}

			convert_to_conservative(params->gamma, p(i, j), vx(i, j), vy(i, j), rho(i, j), m(i, j), impx(i, j), impy(i, j), e(i, j));
			
		}
	}

}
void init_krujok(int local_Nx, int local_Ny, Params* params, vec x, vec y, Field2D<double>& p,
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
	
	stena = 0.3;
	
	
	// ne toro test 1

	p1 = 1.0; // атмосферное давление
	vx1 = 0.75;
	vy1 = 0.0;
	rho1 = 1.0; // плотность воздуха

	p2 = 1.0; // атмосферное давление
	vx2 = 0.;
	vy2 = 0.0;
	rho2 = 0.125; // плотность воздуха
	
	for (int i = 0; i < local_Nx + 2 * params->fict; i++) {
		for (int j = 0; j < local_Ny + 2 * params->fict; j++) {

			d = (x[i] - 0.5) * (x[i] - 0.5) + (y[j] - 0.5) * (y[j] - 0.5);
			
			
			if (d < R*R) {
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
			convert_to_conservative(params->gamma, p(i, j), vx(i, j), vy(i, j), rho(i, j), m(i, j), impx(i, j), impy(i, j), e(i, j));
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

// minmod functions: 0 - Kolgan,1972; 1 - Kolgan,1975; 2 (or another...) - Osher,1984
double minmod(double a, double b, int func) // func - type of minmod we use
{
	if (func == 0)
		return (abs(a) < abs(b) ? a : b);
	else if (func == 1)
	{
		double c = (a + b) / 2;
		return ((abs(a) < abs(b) && abs(a) < abs(c)) ? a : (abs(b) < abs(a) && abs(b) < abs(c) ? b : c));
	}
	else
		return (a * b < 0 ? 0 : (a * a < a * b ? a : b));
}

// is vacuum created?
bool is_vacuum(double gamma, double CL, double UL, double CR, double UR)// C - sound's speed
{
	if (2.0 / (gamma - 1.0) * (CL + CR) < UR - UL) return true;
	return false;
}

// in F calculate U on contact, in FD calculate dU/dP on contact
void prefun(double gamma, double& F, double& FD, double P, double DK, double PK, double CK) // PK  =  PR or PL
{
	if (P <= PK) //rarefaction wave
	{
		F = 2. / (gamma - 1) * CK * (pow((P / PK), (gamma - 1.) / 2. / gamma) - 1.);
		FD = CK * gamma / PK * pow((P / PK), (gamma - 1) / gamma / 2. - 1);
	}
	else //shock wave
	{
		F = (P - PK) / sqrt(DK / 2 * ((gamma + 1) * P + (gamma - 1) * PK));
		FD = 1 / sqrt(2 * DK) * ((gamma + 1) * P + (3 * gamma - 1) * P) / pow(((gamma + 1) * P + (gamma - 1) * PK), 1.5);
	}
}

void guess_p(double gamma, double Quser, double PL, double DL, double UL, double PR, double DR, double UR, double& P)
{
	double CUP, GEL, GER, PMAX, PMIN, PPV, QMAX, EPS = 1.e-8;
	double CL = sqrt(gamma * PL / DL);
	double CR = sqrt(gamma * PR / DR);

	CUP = 0.25 * (DL + DR) * (CL + CR);

	PPV = 0.5 * (PL + PR) + 0.5 * (UL - UR) * CUP; // initial condition from linear task
	PPV = std::max(0.0, PPV);

	PMIN = std::min(PL, PR);
	PMAX = std::max(PL, PR);

	QMAX = PMAX / PMIN;

	double PM, UM;

	if (QMAX <= Quser && (PMIN <= PPV && PPV <= PMAX) || abs(PMIN - PPV) < EPS || abs(PMAX - PPV) < EPS) {
		PM = PPV;
	}
	else {
		if (PPV < PMIN) // two rarefaction waves
		{
			/*

			double PQ = pow(PL / PR, (gamma - 1.0) / (2.0 * gamma));
			UM = (PQ * UL / CL + UR / CR + 2.0 / (gamma - 1.0) * (PQ - 1.0f) / (PQ / CL + 1.0f / CR));
			double PTL = 1.0 + (gamma - 1.0) / 2.0 * (UL - UM) / CL;
			double PTR = 1.0 + (gamma - 1.0) / 2.0 * (UM - UR) / CR;
			PM = 0.5 * (PL * pow(PTL, 2.0 * gamma / (gamma - 1.0)) + PR * pow(PTR, 2.0 * gamma / (gamma - 1.0)));
			*/
			double gamma1 = 0.5 * (gamma - 1.0) / gamma;
			PM = pow(abs(CL + CR - 0.5 * (gamma - 1.0) * (UR - UL)) / abs(CL / pow(PL, gamma1) + CR / pow(PR, gamma1)), 1.0 / gamma1);
		}
		else //two shock waves
		{
			double gamma1 = 2.0 / (gamma + 1.0);
			double gamma2 = (gamma - 1.) / (gamma + 1.);

			GEL = sqrt(gamma1 / DL / (gamma2 * PL + PPV));
			GER = sqrt(gamma1 / DR / (gamma2 * PR + PPV));
			PM = abs(GEL * PL + GER * PR - (UR - UL)) / (GEL + GER);
		}
	}
	P = PM;
}

//  to compute the solution for pressure and velocity on contact
void starpu(double gamma, double Quser, double PL, double DL, double UL, double PR, double DR, double UR, double& UM, double& PM)
{
	int NRITER = 30;
	double CHANGE = 1e6, FL, FLD, FR, FRD, POLD, TOLPRE = 1.0e-6, UDIFF;
	double CL = sqrt(gamma * PL / DL);
	double CR = sqrt(gamma * PR / DR);


	if (is_vacuum(gamma, CL, UL, CR, UR)) {
		std::cout << "Vacuum is generated" << std::endl << "Program stopped" << std::endl;
		exit(228);
	}
	guess_p(gamma, Quser, PL, DL, UL, PR, DR, UR, POLD);		// initial condition for pressure

	if (POLD < 0.)
	{
		std::cout << "Pressure is negative" << std::endl;
		exit(105);
	}

	if (POLD != POLD)
	{
		std::cout << "Pressure is NaN after guess_p" << std::endl;
		exit(-1488);
	}

	UDIFF = UR - UL;

	for (int i = 0; i < NRITER; ++i)
	{
		prefun(gamma, FL, FLD, POLD, DL, PL, CL);
		prefun(gamma, FR, FRD, POLD, DR, PR, CR);

		if (POLD < 0.)
		{
			std::cout << "Pressure is negative" << std::endl;
			exit(106);
		}
		if (POLD != POLD)
		{
			std::cout << "Pressure is NaN in iter" << std::endl;
			exit(1323415);
		}

		PM = POLD - (FL + FR + UDIFF) / (FLD + FRD);
		CHANGE = 2.0 * abs((PM - POLD) / (PM + POLD));
		if (CHANGE < TOLPRE) break;
	}
	UM = 0.5 * (UL + UR + FR - FL);

	return;
}

// to sample the solution throughout the wave pattern, PM, UM - contact's parameters
void sample(double gamma, double PL, double DL, double UL, double PR, double DR, double UR,
	double UM, double PM, double S, double& D, double& U, double& P)
{
	double C, CML, CMR, PML, PMR, SHL, SHR, SL, SR, STL, STR;
	//
	// C - local sound speed in rarefaction wave
	// CML, CMR - sound speed on the left / right side of contact
	// SHL, SHR - speeds of head of left/right rarefaction wave
	// SL, SR - left/right shock speed
	// STL, STR - speeds of tail of left/right rarefaction wave
	//
	double CL = sqrt(gamma * PL / DL);
	double CR = sqrt(gamma * PR / DR);

	if (S <= UM) // point is left from contact
	{
		if (PM <= PL) //left rarefaction wave
		{
			SHL = UL - CL;
			if (S <= SHL) //point is left from rarefaction wave - left parameters
			{
				D = DL;
				U = UL;
				P = PL;
			}
			else
			{
				CML = CL * pow((PM / PL), (gamma - 1.0) / (2.0 * gamma));
				STL = UM - CML;
				if (S > STL) // point is a state with *
				{
					D = DL * pow((PM / PL), (1. / gamma));
					U = UM;
					P = PM;
				}
				else //point is in a left rarefaction wave
				{
					U = 2.0 / (gamma + 1.0) * (CL + (gamma - 1.0) / 2.0 * UL + S);
					C = 2.0 / (gamma + 1.0) * (CL + (gamma - 1.0) / 2.0 * (UL - S));
					D = DL * pow((C / CL), 2.0 / (gamma - 1.0));
					P = PL * pow((C / CL), 2.0 * gamma / (gamma - 1.0));
				}
			}
		}
		else // left shock wave
		{
			PML = PM / PL;
			SL = UL - CL * sqrt((gamma + 1.0) / (2.0 * gamma) * PML + (gamma - 1.0) / (2.0 * gamma));
			if (S <= SL) //point is left from shock wave - left parameters
			{
				D = DL;
				U = UL;
				P = PL;
			}
			else // point is a state with *
			{
				D = DL * (PML + (gamma - 1.0) / (gamma + 1.0)) / (PML * (gamma - 1.0) / (gamma + 1.0) + 1.);
				U = UM;
				P = PM;
			}
		}
	}
	else // point is right from contact
	{
		if (PM > PR) //right from shock wave
		{
			PMR = PM / PR;
			SR = UR + CR * sqrt((gamma + 1.0) / (2.0 * gamma) * PMR + (gamma - 1.0) / (2.0 * gamma));
			if (S >= SR) //point is right from shock wave - right parameters
			{
				D = DR;
				U = UR;
				P = PR;
			}
			else // point is a state with *
			{
				D = DR * (PMR + (gamma - 1.0) / (gamma + 1.0)) / (PMR * (gamma - 1.0) / (gamma + 1.0) + 1.);
				U = UM;
				P = PM;
			}
		}
		else // right rarefaction wave
		{
			SHR = UR + CR;
			if (S >= SHR) //point is right from rarefaction wave - right parameters
			{
				D = DR;
				U = UR;
				P = PR;
			}
			else
			{
				CMR = CR * pow((PM / PR), (gamma - 1.0) / (2.0 * gamma));
				STR = UM + CMR;
				if (S <= STR) // point is a state with *
				{
					D = DR * pow((PM / PR), (1. / ((gamma - 1.0) / (2.0 * gamma))));
					U = UM;
					P = PM;
				}
				else // point is in a right rarefaction wave
				{
					U = 2.0 / (gamma + 1.0) * (-CR + (gamma - 1.0) / 2.0 * UR + S);
					C = 2.0 / (gamma + 1.0) * (CR - (gamma - 1.0) / 2.0 * (UR - S));
					D = DR * pow((C / CR), 2.0 / (gamma - 1.0));
					P = PR * pow((C / CR), 2.0 * gamma / (gamma - 1.0));
				}
			}
		}
	}
	return;
}


// Riemann solver
void Riemann_solver(double gamma, double Quser, double PL, double DL, double UL, double PR, double DR, double UR,
	double& D, double& U, double& P)
{
	double CL = sqrt(gamma * PL / DL); // left sound speed
	double CR = sqrt(gamma * PR / DR); // right sound speed
	double PM, UM; // pressure and velocity on contact
	double S = 0.;

	// Проверка входных данных
	if (DL <= 0.0 || DR <= 0.0 || PL <= 0.0 || PR <= 0.0) {
		std::cerr << "Invalid input to Riemann solver" << std::endl;
		exit(1);
	}

	// vacuum generation
	if (is_vacuum(gamma, CL, UL, CR, UR)) {
		std::cout << "Vacuum is generated" << std::endl << "Program stopped" << std::endl;
		exit(228);
	}
	// iteration pressure and velocity
	//void starpu(double gamma, double Quser, double PL, double DL, double UL, double PR, double DR, double UR, double& UM, double& PM)

	starpu(gamma, Quser, PL, DL, UL, PR, DR, UR, UM, PM);

	// found results
	sample(gamma, PL, DL, UL, PR, DR, UR, UM, PM, 0., D, U, P);

	if (P <= 0.0 || D <= 0.0 || !std::isfinite(P) || !std::isfinite(U)) {
		std::cerr << "Invalid output from Riemann solver" << std::endl;
		exit(1);
	}
}

void Godunov_flux_x(double gamma, double p, double vx, double vy, double rho, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
	double m, impx, impy, e;
	
	convert_to_conservative(gamma, p, vx, vy, rho, m, impx, impy, e);
	Fm = rho * vx;
	Fimpx = Fm * vx + p;
	Fimpy = Fm * vy;
	Fe = (p + e) * vx;
}

void Godunov_flux_y(double gamma, double p, double vx, double vy, double rho, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
	double m, impx, impy, e;

	convert_to_conservative(gamma, p, vx, vy, rho, m, impx, impy, e);
	Fm = rho * vy;
	Fimpx = Fm * vx;
	Fimpy = Fm * vy + p;
	Fe = (p + e) * vy;
}

void Godunov_method_x(double gamma, double Quser, double ml, double impxl, double impyl, double el,
	double mr, double impxr, double impyr, double er,
	double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
	double p, vx, vy, rho;
	double pl, vxl, vyl, rhol;
	double pr, vxr, vyr, rhor;

	convert_from_conservative(gamma, pl, vxl, vyl, rhol, ml, impxl, impyl, el);
	convert_from_conservative(gamma, pr, vxr, vyr, rhor, mr, impxr, impyr, er);

	Riemann_solver(gamma, Quser, pl, rhol, vxl, pr, rhor, vxr, rho, vx, p);

	if (vx >= 0)
		vy = impyl / ml;
	else
		vy = impyr / mr;

	Godunov_flux_x(gamma, p, vx, vy, rho, Fm, Fimpx, Fimpy, Fe);
}

void Godunov_method_y(double gamma, double Quser, double md, double impxd, double impyd, double ed,
	double mu, double impxu, double impyu, double eu,
	double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
	double p, vx, vy, rho;
	double pd, vxd, vyd, rhod;
	double pu, vxu, vyu, rhou;

	convert_from_conservative(gamma, pd, vxd, vyd, rhod, md, impxd, impyd, ed);
	convert_from_conservative(gamma, pu, vxu, vyu, rhou, mu, impxu, impyu, eu);

	Riemann_solver(gamma, Quser, pd, rhod, vyd, pu, rhou, vyu, rho, vy, p);

	if (vy >= 0)
		vx = impxd / md;
	else
		vx = impxu / mu;

	Godunov_flux_y(gamma, p, vx, vy, rho, Fm, Fimpx, Fimpy, Fe);
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

			e_next(i, j) +=  (wx + wy) / (params->gamma - 1.0);
		}
	}

}
void GK_2d(int argc, char** argv)
{
	int rank, size;
	MPI_Init(&argc, &argv);
	double start = MPI_Wtime();
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
	
	Params* params = new Params(std::stol(argv[1]), std::stol(argv[2]), 1, 0.15, 4.0, 2.0, std::stod(argv[4]), 1.4, 2.0);
	
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

	vec xc(local_Nx + 2 * params->fict);
	vec x(local_Nx + 1 + 2 * params->fict);
	vec yc(local_Ny + 2 * params->fict);
	vec y(local_Ny + 1 + 2 * params->fict);

	//+ 2 для оптимизации обмена между процессами
	Field2D<double> p(local_Nx + 2 * params->fict, local_Ny + 2 * params->fict), vx(local_Nx + 2 * params->fict, local_Ny + 2 * params->fict), vy(local_Nx + 2 * params->fict, local_Ny + 2 * params->fict), rho(local_Nx + 2 * params->fict, local_Ny + 2 * params->fict);
	Field2D<double> m(local_Nx + 2 * params->fict, local_Ny + 2 * params->fict), impx(local_Nx + 2 * params->fict, local_Ny + 2 * params->fict), impy(local_Nx + 2 * params->fict, local_Ny + 2 * params->fict), e(local_Nx + 2 * params->fict, local_Ny + 2 * params->fict);
	Field2D<double> m_next(local_Nx + 2 * params->fict, local_Ny + 2 * params->fict), impx_next(local_Nx + 2 * params->fict, local_Ny + 2 * params->fict), impy_next(local_Nx + 2 * params->fict, local_Ny + 2 * params->fict), e_next(local_Nx + 2 * params->fict, local_Ny + 2 * params->fict);

	// задание границ клина
	Field2D<int> klin(local_Nx + 2 * params->fict, local_Ny + 2 * params->fict);

	double mb, impxb, impyb, eb, pb, vxb, vyb, rhob;
	double FmL, FimpxL, FimpyL, FeL, FmR, FimpxR, FimpyR, FeR;

	double pl, vl, rhol, pr, vr, rhor;
	double ml, impxl, impyl, el, mr, impxr, impyr, er;
	double dm, dimp, de;

	int step = 0, max_step = 1000000;// 201;
	double time = 0.0;
	double dt_temp;

	int out = 10;
	int minmod_type = 1, bound_left = 2, bound_right = 1, bound_down = 0, bound_up = 0;
	int viscosity_type = 2;

	int left_b = 0, down_b = 0;

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

		
	grid(params, x, y, xc, yc, local_Nx, local_Ny, left_b * (params->dx), down_b * (params->dy));
	
	init_mach(local_Nx, local_Ny, params, x, y, p, vx, vy, rho, m, impx, impy, e, klin);
	update_boundaries(grid_comm, x, y, p, vx, vy, rho, m, impx, impy, e, klin, local_Nx, local_Ny,
			left, right, up, down,
			bound_left, bound_right, bound_down, bound_up,
			params);
		
	data_to_file(local_Nx, local_Ny, params->fict, params->gamma, 0.0, x, y, p, vx, vy, rho, rank, dims, (argv[3]));
	
	
	while (time < params->T && step < max_step)
	{	
		
		update_boundaries(grid_comm, x, y, p, vx, vy, rho, m, impx, impy, e, klin, local_Nx, local_Ny,
			left, right, up, down,
			bound_left, bound_right, bound_down, bound_up,
			params);

		
		dt_temp = get_dt(params, p, rho, vx, vy, x, y, local_Nx, local_Ny);
		
		double dt = 0.0;

		MPI_Allreduce(&dt_temp, &dt, 1, MPI_DOUBLE, MPI_MIN, grid_comm);

		time += dt;
		
		/*   START CALC   */
		for (int i = params->fict; i < local_Nx + (params->fict); i++)
		{
			for (int j = params->fict; j < local_Ny + (params->fict); j++)
			{
				//		расчет вдоль х
				
				//		поток через левую грань ячейки	

				{
					ml = m(i - 1, j) + 0.5 * minmod(m(i - 1, j) - m(i - 2, j), m(i, j) - m(i - 1, j), minmod_type);
					impxl = impx(i - 1, j) + 0.5 * minmod(impx(i - 1, j) - impx(i - 2, j), impx(i, j) - impx(i - 1, j), minmod_type);
					impyl = impy(i - 1, j) + 0.5 * minmod(impy(i - 1, j) - impy(i - 2, j), impy(i, j) - impy(i - 1, j), minmod_type);
					el = e(i - 1, j) + 0.5 * minmod(e(i - 1, j) - e(i - 2, j), e(i, j) - e(i - 1, j), minmod_type);

					mr = m(i, j) - 0.5 * minmod(m(i, j) - m(i - 1, j), m(i + 1, j) - m(i, j), minmod_type);
					impxr = impx(i, j) - 0.5 * minmod(impx(i, j) - impx(i - 1, j), impx(i + 1, j) - impx(i, j), minmod_type);
					impyr = impy(i, j) - 0.5 * minmod(impy(i, j) - impy(i - 1, j), impy(i + 1, j) - impy(i, j), minmod_type);
					er = e(i, j) - 0.5 * minmod(e(i, j) - e(i - 1, j), e(i + 1, j) - e(i, j), minmod_type);
				}
				
				Godunov_method_x(params->gamma, params->Quser, ml, impxl, impyl, el, mr, impxr, impyr, er, FmL, FimpxL, FimpyL, FeL);
				
				//		поток через правую грань ячейки
				
				{
					ml = m(i, j) + 0.5 * minmod(m(i, j) - m(i - 1, j), m(i + 1, j) - m(i, j), minmod_type);
					impxl = impx(i, j) + 0.5 * minmod(impx(i, j) - impx(i - 1, j), impx(i + 1, j) - impx(i, j), minmod_type);
					impyl = impy(i, j) + 0.5 * minmod(impy(i, j) - impy(i - 1, j), impy(i + 1, j) - impy(i, j), minmod_type);
					el = e(i, j) + 0.5 * minmod(e(i, j) - e(i - 1, j), e(i + 1, j) - e(i, j), minmod_type);

					mr = m(i + 1, j) - 0.5 * minmod(m(i + 1, j) - m(i, j), m(i + 2, j) - m(i + 1, j), minmod_type);
					impxr = impx(i + 1, j) - 0.5 * minmod(impx(i + 1, j) - impx(i, j), impx(i + 2, j) - impx(i + 1, j), minmod_type);
					impyr = impy(i + 1, j) - 0.5 * minmod(impy(i + 1, j) - impy(i, j), impy(i + 2, j) - impy(i + 1, j), minmod_type);
					er = e(i + 1, j) - 0.5 * minmod(e(i + 1, j) - e(i, j), e(i + 2, j) - e(i + 1, j), minmod_type);
				}
				
				Godunov_method_x(params->gamma, params->Quser, ml, impxl, impyl, el, mr, impxr, impyr, er, FmR, FimpxR, FimpyR, FeR);

				m_next(i, j) = m(i, j) - dt * (FmR - FmL) / params->dx;
				impx_next(i, j) = impx(i, j) - dt * (FimpxR - FimpxL) / params->dx;
				impy_next(i, j) = impy(i, j) - dt * (FimpyR - FimpyL) / params->dx;
				e_next(i, j) = e(i, j) - dt * (FeR - FeL) / params->dx;

				//		расчёт вдоль y
				
				//		поток через нижнюю грань
				
				
				{
					ml = m(i, j - 1) + 0.5 * minmod(m(i, j - 1) - m(i, j - 2), m(i, j) - m(i, j - 1), minmod_type);
					impxl = impx(i, j - 1) + 0.5 * minmod(impx(i, j - 1) - impx(i, j - 2), impx(i, j) - impx(i, j - 1), minmod_type);
					impyl = impy(i, j - 1) + 0.5 * minmod(impy(i, j - 1) - impy(i, j - 2), impy(i, j) - impy(i, j - 1), minmod_type);
					el = e(i, j - 1) + 0.5 * minmod(e(i, j - 1) - e(i, j - 2), e(i, j) - e(i, j - 1), minmod_type);

					mr = m(i, j) - 0.5 * minmod(m(i, j) - m(i, j - 1), m(i, j + 1) - m(i, j), minmod_type);
					impxr = impx(i, j) - 0.5 * minmod(impx(i, j) - impx(i, j - 1), impx(i, j + 1) - impx(i, j), minmod_type);
					impyr = impy(i, j) - 0.5 * minmod(impy(i, j) - impy(i, j - 1), impy(i, j + 1) - impy(i, j), minmod_type);
					er = e(i, j) - 0.5 * minmod(e(i, j) - e(i, j - 1), e(i, j + 1) - e(i, j), minmod_type);
				}
				
				Godunov_method_y(params->gamma, params->Quser, ml, impxl, impyl, el, mr, impxr, impyr, er, FmL, FimpxL, FimpyL, FeL);
				
				//		поток через верхнюю грань ячейки
				
				{
					ml = m(i, j) + 0.5 * minmod(m(i, j) - m(i, j - 1), m(i, j + 1) - m(i, j), minmod_type);
					impxl = impx(i, j) + 0.5 * minmod(impx(i, j) - impx(i, j - 1), impx(i, j + 1) - impx(i, j), minmod_type);
					impyl = impy(i, j) + 0.5 * minmod(impy(i, j) - impy(i, j - 1), impy(i, j + 1) - impy(i, j), minmod_type);
					el = e(i, j) + 0.5 * minmod(e(i, j) - e(i, j - 1), e(i, j + 1) - e(i, j), minmod_type);

					mr = m(i, j + 1) - 0.5 * minmod(m(i, j + 1) - m(i, j), m(i, j + 2) - m(i, j + 1), minmod_type);
					impxr = impx(i, j + 1) - 0.5 * minmod(impx(i, j + 1) - impx(i, j), impx(i, j + 2) - impx(i, j + 1), minmod_type);
					impyr = impy(i, j + 1) - 0.5 * minmod(impy(i, j + 1) - impy(i, j), impy(i, j + 2) - impy(i, j + 1), minmod_type);
					er = e(i, j + 1) - 0.5 * minmod(e(i, j + 1) - e(i, j), e(i, j + 2) - e(i, j + 1), minmod_type);
				}
				
				Godunov_method_y(params->gamma, params->Quser, ml, impxl, impyl, el, mr, impxr, impyr, er, FmR, FimpxR, FimpyR, FeR);
			
				
				m_next(i, j) = m_next(i, j) - dt * (FmR - FmL) / params->dy;
				impx_next(i, j) = impx_next(i, j) - dt * (FimpxR - FimpxL) / params->dy;
				impy_next(i, j) = impy_next(i, j) - dt * (FimpyR - FimpyL) / params->dy;
				e_next(i, j) = e_next(i, j) - dt * (FeR - FeL) / params->dy;
				
			}
		}
		
		update_boundaries(grid_comm, x, y, p, vx, vy, rho, m, impx, impy, e, klin, local_Nx, local_Ny,
							left, right, up, down, bound_left, bound_right,
							bound_down, bound_up, params);

		if (viscosity_type != 0)
		{
			Artificial_viscosity(viscosity_type, p, vx, vy, rho, m_next, impx_next, impy_next, e_next, local_Nx, local_Ny, params);
		}

		/*    END CALC    */

		// update parameters
		for (int i = params->fict; i < local_Nx + params->fict; ++i)
		{
			for (int j = params->fict; j < local_Ny + params->fict; ++j)
			{
				
				if (klin(i, j) == 0) {
				
					m(i, j) = m_next(i, j);
					impx(i, j) = impx_next(i, j);
					impy(i, j) = impy_next(i, j);
					e(i, j) = e_next(i, j);
					
					convert_from_conservative(params->gamma, p(i, j), vx(i, j), vy(i, j), rho(i, j), m(i, j), impx(i, j), impy(i, j), e(i, j));
				}
			}
		}

		if (step % out == 0)
			data_to_file(local_Nx, local_Ny, params->fict, params->gamma, time, x, y, p, vx, vy, rho, rank, dims, (argv[3]));
		step += 1;
		//if (rank == 0 && !(step % out)) { 	std::cout << step << std::endl;}

	}
	double end = MPI_Wtime();
	if (rank == 0) {
	
	printf("proc,time,N,M\n");
	printf("%d,%f,%d,%d\n", rank,end - start, params->Nx, params->Ny);
	
	}
	else {

	printf("%d,%f,%d,%d\n", rank, end - start, params->Nx, params->Ny);
	}
	MPI_Finalize();
	return;
}


int main(int argc, char** argv)
{
	GK_2d(argc, argv);

	return 0;
}