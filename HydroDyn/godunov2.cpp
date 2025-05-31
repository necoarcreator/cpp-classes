# include "godunov2.h"
#include <numeric>
godunov::godunov(string file_path): krootoy(file_path) {
	imp = new Grid(N_x);
	imp_new = new Grid(N_x);
    p_an = new Grid(N_x);
    v_an = new Grid(N_x);
    rho_an = new Grid(N_x);
    m_an = new Grid(N_x);
    imp_an = new Grid(N_x);
    e_an = new Grid(N_x);

	delete x;
	
	x = new Grid(N_x);
	m_new = new Grid(N_x); // Масса

};


godunov::~godunov()  {
	
	delete imp;
	
	delete imp_new;

    delete p_an;

    delete v_an;

    delete rho_an;

    delete m_an;

    delete imp_an;

    delete e_an;
};

// граничные условия унаследованы

// Get_dt (), makeGrid () унаследованы

// writeResults унаследованы

// pokushai тесты унаследованы

void godunov::Convert_cons_to_noncons(double& p, double& v, double& rho, double m, double imp, double e) {
    p = (gamma - 1.0) * (e - 0.5 * (pow(imp, 2.0)) / m);
    v = imp / m;
    rho = m;
}

void godunov::Convert_cons_to_noncons(Grid* p, Grid* v, Grid* rho, Grid* m, Grid* imp, Grid* e) {
    for (int i = 0; i < N_x; i++) {
        (*p)[i] = (gamma - 1.0) * ((*e)[i] - 0.5 * (pow((*imp)[i], 2.0)) / (*m)[i]);
        (*v)[i] = (*imp)[i] / (*m)[i];
        (*rho)[i] = (*m)[i];
    }
}

void godunov::Convert_noncons_to_cons(double p, double v, double rho, double& m, double& imp, double& e) {
    m = rho;
    imp = rho * v;
    e = 0.5 * rho * (pow(v, 2.0)) + p / (gamma - 1.0);
}

void godunov::Convert_noncons_to_cons(Grid* p, Grid* v, Grid* rho, Grid* m, Grid* imp, Grid* e) {
    for (int i = 0; i < N_x; i++) {
        (*m)[i] = (*rho)[i];
        (*imp)[i] = (*rho)[i] * (*v)[i];
        (*e)[i] = 0.5 * (*rho)[i] * (pow((*v)[i], 2.0)) + (*p)[i] / (gamma - 1.0);
    }
}

double godunov::Get_dt () {
	double new_step = 10000000;
	double c; // скорость звука
	double t_step; // текущий шаг по времени
	double p_dt = 0., rho_dt = 0., v_dt = 0.;
	
	
	for (int i = fict; i < N_x - fict; ++i) {
		Convert_cons_to_noncons (p_dt, v_dt, rho_dt, (*m)[i], (*imp)[i], (*e)[i]);
		c = Sound_velocity (p_dt, rho_dt);
		t_step = CFL * ((*xb)[i+1]-(*xb)[i])/(abs(v_dt)+c);
		
		if (t_step < new_step) new_step = t_step; 
	}
	
	return new_step;
}

double godunov::Sound_velocity (double p_sound_velocity, double rho_sound_velocity) {
	return sqrt (gamma * p_sound_velocity/rho_sound_velocity); }

void godunov::Solver_godunov(int test, string file_path_res, string file_path_an, string format, char _delim, int fo) {
	double imp_r = 0.;  // импульс справа
	double imp_l = 0.; // импульс слева
	double p_r = 0.; // давление справа
	double p_l = 0.; // давление слева
	double v_r = 0.; // скорость справа
	double v_l = 0.; // скорость слева
	double m_l = 0., e_l = 0., m_r = 0., e_r = 0.; // масса и полныая энергия
	double rho_l = 0., rho_r = 0.; //плотности слева и справа
	double dt, dx;
	
	double F_m_l = 0., F_imp_l = 0., F_e_l = 0., F_m_r = 0., F_imp_r = 0., F_e_r = 0.; // потоки консервативные
	
	double time = 0;
	
	int step_number = 0;
	
	makeGrid();

    Grid* p_an_start = new Grid(N_x);
    Grid* v_an_start = new Grid(N_x);
    Grid* rho_an_start = new Grid(N_x);
    Grid* m_an_start = new Grid(N_x);
    Grid* imp_an_start = new Grid(N_x);
    Grid* e_an_start = new Grid(N_x);

    fict = 0;

	switch(test) {
        case 0:
            Test0(N_x, p, rho, v);
            Test0(N_x, p_an_start, rho_an_start, v_an_start);
            break;
        case 1:
            Test1(N_x, p, rho, v);
            Test1(N_x, p_an_start, rho_an_start, v_an_start);
            break;
        case 2:
            Test2(N_x, p, rho, v);
            Test2(N_x, p_an_start, rho_an_start, v_an_start);
            break;
        case 3:
            Test3(N_x, p, rho, v);
            Test3(N_x, p_an_start, rho_an_start, v_an_start);
            break;
        case 4:
            Test4(N_x, p, rho, v);
            Test4(N_x, p_an_start, rho_an_start, v_an_start);
            break;
        case 5:
            Test5(N_x, p, rho, v);
            Test5(N_x, p_an_start, rho_an_start, v_an_start);
            break;
        default:
            std::cout << "Invalid test number" << std::endl;
            exit(1);
            break;
    }
	
    Convert_noncons_to_cons(p, v, rho, m, imp, e);
    Convert_noncons_to_cons(p_an_start, v_an_start, rho_an_start, m_an_start, imp_an_start, e_an_start);
	
	writeResults(file_path_res, file_path_an, format, _delim, time, test);
	
    double end_time = T_1 - T_0;

	while (time < end_time) {
		dt = Get_dt();
		for (int i = 0; i < N_x; ++i) {
			if (i != 0) {
				m_l = (*m)[i-1];
				imp_l = (*imp) [i-1];
				e_l = (*e) [i-1];
				m_r = (*m)[i];
				imp_r = (*imp) [i];
				e_r = (*e) [i];
			} else {
                BoundaryOpenCons(0, *m, *imp, *e);
                m_l = (*m)[0];
                imp_l = (*imp)[0];
                e_l = (*e)[0];
				m_r = (*m)[i];
				imp_r = (*imp) [i];
				e_r = (*e) [i];
			}
			
			Godunov_flux(m_l, imp_l, e_l, m_r, imp_r, e_r, F_m_l, F_imp_l, F_e_l);
			
			if (i != N_x-1) {
				m_l = (*m)[i];
				imp_l = (*imp) [i];
				e_l = (*e) [i];
				m_r = (*m)[i+1];
				imp_r = (*imp) [i+1];
				e_r = (*e) [i+1];	
			} else {
				BoundaryOpenCons(1, *m, *imp, *e);
                m_r = (*m)[N_x - 1];
                imp_r = (*imp)[N_x - 1];
                e_r = (*e)[N_x - 1];
				m_l = (*m)[i];
				imp_l = (*imp) [i];
				e_l = (*e) [i];
            }
			Godunov_flux(m_l, imp_l, e_l, m_r, imp_r, e_r, F_m_r, F_imp_r, F_e_r);
			
			dx = (Lx_1 - Lx_0)/(N_x);
			
			(*m_new)[i] = (*m)[i] - dt * (F_m_r - F_m_l)/dx;
			(*imp_new)[i] = (*imp)[i] - dt * (F_imp_r - F_imp_l)/dx;
			(*e_new)[i] = (*e)[i] - dt * (F_e_r - F_e_l)/dx;
		}
			
		for (int i = 0; i < N_x; ++i) {
			(*imp)[i] = (*imp_new)[i];
			(*m)[i] = (*m_new)[i];
			(*e)[i] = (*e_new)[i];
			Convert_cons_to_noncons ((*p)[i], (*v)[i], (*rho)[i], (*m)[i], (*imp)[i], (*e)[i]);
		}
		
		time += dt;
		step_number++;
        if (step_number % fo == 0) {
            double l_coord = (*x)[N_x / 2];

            for (size_t i = 0; i < N_x - 1; i++) {
                Analytic(time, l_coord, (*x)[i], (*m_an_start)[0], (*imp_an_start)[0], (*e_an_start)[0], (*m_an_start)[N_x - 1], 
                     (*imp_an_start)[N_x - 1], (*e_an_start)[N_x - 1], (*p_an)[i], (*v_an)[i], (*rho_an)[i], (*m_an)[i], (*imp_an)[i], (*e_an)[i]);
            }
            writeResults(file_path_res, file_path_an, format, _delim, time, test);
        }
	}
    delete p_an_start;
    delete v_an_start;
    delete rho_an_start;
    delete m_an_start;
    delete imp_an_start;
    delete e_an_start;
}

void godunov::Godunov_flux(double m_l, double imp_l, double e_l, double m_r, double imp_r, double e_r, double& Fm, double& Fimp_x, double& Fe) {
	double p_flux = 0., v_flux = 0., rho_flux = 0.;
    double pl = 0., vl = 0., rhol = 0.;
    double pr = 0., vr = 0., rhor = 0.;
	
    Convert_cons_to_noncons(pl, vl, rhol, m_l, imp_l, e_l);
	Convert_cons_to_noncons(pr, vr, rhor, m_r, imp_r, e_r);
	
    /* расчет потока Годунова по вектору неконсервативных переменных*/
    Diff_flux_ncons_x(p_flux, v_flux, rho_flux, Fm, Fimp_x, Fe);
}

void godunov::Diff_flux_ncons_x(double p_flux, double v_flux, double rho, double& Fm, double& Fimp_x, double& Fe) {
    double m_flux = 0., impx = 0., e_flux = 0.; /* консервативные переменные */

    Convert_noncons_to_cons(p_flux, v_flux, rho, m_flux, impx, e_flux);
	
    Fm = rho * v_flux;
    Fimp_x = Fm * v_flux + p_flux;
    Fe = (p_flux + e_flux) * v_flux;
}

void godunov::Riman_solver(double rhol, double vl, double pl, double rhor, double vr, double pr, double& p_res, double& v_res, double& r_res) {
    double cr, cl;
    double p_cont, v_cont; //значения, полученные итерационной сшвкой на контактном разрыве
    double p_riman, v_riman, r_riman; //значения, воссстановленные на границе по решению задачи Римана

    cl = Sound_velocity(pl, rhol);
    cr = Sound_velocity(pr, rhor);

    if (2.0 * (cl + cr) / (gamma - 1.0) <= vr - vl) {
        /* случай возникновения вакуума */
        printf("\nContact_pressure_velocity -> vacuum is generated\n");
    }
    /* итерационная процедура расчета давления и скорости газа на контактном разрыве*/
    Contact_pressure_velocity(pl, vl, rhol, cl, pr, vr, rhor, cr, p_cont, v_cont);
    /* отбор решения */
    Sample_solid_solution(pl, vl, rhol, cl, pr, vr, rhor, cr, p_cont, v_cont, 0.0, p_riman, v_riman, r_riman);
    p_res = p_riman;
    v_res = v_riman;
    r_res = r_riman;
}

/* отбор решения */
void godunov::Sample_solid_solution(double pl, double vl, double rhol, double cl, double pr, double vr, double rhor, double cr,
    double p_cont, double v_cont, double s, double& p_res, double& v_res, double& rho_res) {

    double g1, g2, g3, g4, g5, g6, g7;      /* вспомогательные переменные, производные от показателя адиабаты */

                                               /* скорости левых волн */
    double shl, stl;        /* скорости "головы" и "хвоста" левой волны разрежения */
    double sl;              /* скорость левой ударной волны */

    /* скорости правых волн */
    double shr, str;        /* скорости "головы" и "хвоста" правой волны разрежения */
    double sr;              /* скорость правой ударной волны */

    double cml, cmr;        /* скорости звука слева и справа от контактного разрыва */
    double c;               /* локальная скорость звука внутри волны разрежения */
    double p_ratio;
    double r_solution, v_solution, p_solution;         /* отобранные значения объемной доли, плотности, скорости и давления */

    /* производные от показателя адиабаты */
    g1 = 0.5 * (gamma - 1.0) / gamma;
    g2 = 0.5 * (gamma + 1.0) / gamma;
    g3 = 2.0 * gamma / (gamma - 1.0);
    g4 = 2.0 / (gamma - 1.0);
    g5 = 2.0 / (gamma + 1.0);
    g6 = (gamma - 1.0) / (gamma + 1.0);
    g7 = 0.5 * (gamma - 1.0);

    if (s <= v_cont) {
        /* рассматриваемая точка - слева от контактного разрыва */
        if (p_cont <= pl) {
            /* левая волна разрежения */
            shl = vl - cl;
            if (s <= shl) {
                /* параметры слева от разрыва */
                r_solution = rhol;
                v_solution = vl;
                p_solution = pl;
            }
            else {
                cml = cl * pow(p_cont / pl, g1);
                stl = v_cont - cml;
                if (s > stl) {
                    /* параметры слева от контактного разрыва */
                    r_solution = rhol * pow(p_cont / pl, 1.0 / gamma);
                    v_solution = v_cont;
                    p_solution = p_cont;
                }
                else {
                    /* параметры внутри левой волны разрежения */
                    v_solution = g5 * (cl + g7 * vl + s);
                    c = g5 * (cl + g7 * (vl - s));
                    r_solution = rhol * pow(c / cl, g4);
                    p_solution = pl * pow(c / cl, g3);
                }
            }
        }
        else {
            /* левая ударная волна */
            p_ratio = p_cont / pl;
            sl = vl - cl * std::sqrt(g2 * p_ratio + g1);
            if (s <= sl) {
                /* параметры слева от разрыва */
                r_solution = rhol;
                v_solution = vl;
                p_solution = pl;
            }
            else {
                /* параметры за левой ударной волной */
                r_solution = rhol * (p_ratio + g6) / (p_ratio * g6 + 1.0);
                v_solution = v_cont;
                p_solution = p_cont;
            }
        }
    }
    else {
        /* рассматриваемая точка - справа от контактного разрыва */
        if (p_cont > pr) {
            /* правая ударная волна */
            p_ratio = p_cont / pr;
            sr = vr + cr * std::sqrt(g2 * p_ratio + g1);
            if (s >= sr) {
                /* параметры справа от разрыва */
                r_solution = rhor;
                v_solution = vr;
                p_solution = pr;
            }
            else {
                /* параметры за правой ударной волной */
                r_solution = rhor * (p_ratio + g6) / (p_ratio * g6 + 1.0);
                v_solution = v_cont;
                p_solution = p_cont;
            }
        }
        else {
            /* правая волна разрежения */
            shr = vr + cr;
            if (s >= shr) {
                /* параметры справа от разрыва */
                r_solution = rhor;
                v_solution = vr;
                p_solution = pr;
            }
            else {
                cmr = cr * pow(p_cont / pr, g1);
                str = v_cont + cmr;
                if (s <= str) {
                    /* параметры справа от контактного разрыва */
                    r_solution = rhor * pow(p_cont / pr, 1.0 / gamma);
                    v_solution = v_cont;
                    p_solution = p_cont;
                }
                else {
                    /* параметры внутри правой волны разрежения */
                    v_solution = g5 * (-cr + g7 * vr + s);
                    c = g5 * (cr - g7 * (vr - s));
                    r_solution = rhor * pow(c / cr, g4);
                    p_solution = pr * pow(c / cr, g3);
                }
            }
        }
    }
    /* формирование выходного вектора с результатом */
    rho_res = r_solution;
    v_res = v_solution;
    p_res = p_solution;
}

/* Итерационная процедура расчета давления и скорости на контактном разрыве */
void godunov::Contact_pressure_velocity(double pl, double vl, double rhol, double cl,
    double pr, double vr, double rhor, double cr, double& p_cont, double& v_cont) {
    double p_old;         /* значение давления на предыдущей итерации */
    double fl, fr;        /* значения функций */
    double fld, frd;      /* значения производных */
    int iter_num = 0;     /* количество проведенных итераций */
    int iter_max = 300;   /* максимальное количество итераций */
    double criteria;      /* переменная для определения сходимости */
    double eps = 1.e-8;
    if (2.0 * (cl + cr) / (gamma - 1.0) <= vr - vl) {
        /* случай возникновения вакуума */
        printf("\nContact_pressure_velocity -> vacuum is generated\n");
    }
    /* расчет начального приближения для давления */
    p_old = Pressure_initial_guess(pl, vl, rhol, cl, pr, vr, rhor, cr);
    if (p_old < 0.0) {
        printf("\nContact_pressure_velocity -> initial pressure guess is negative ");
    }
    /* решение нелинейного уравнения для нахождения давления на контактном разрыве методом Ньютона-Рафсона */
    do {
        Calc_F_and_DF(p_old, pl, vl, rhol, cl, &fl, &fld);
        Calc_F_and_DF(p_old, pr, vr, rhor, cr, &fr, &frd);
        p_cont = p_old - (fl + fr + vr - vl) / (fld + frd);
        criteria = 2.0 * std::fabs((p_cont - p_old) / (p_cont + p_old));
        iter_num++;
        if (iter_num > iter_max) {
            printf("\nContact_pressure_velocity -> number of iterations exceeds the maximum value.\n");
        }
        if (p_cont < 0.0) {
            printf("\nContact_pressure_velocity -> pressure is negative.\n");
        }
        p_old = p_cont;
    } while (criteria > eps);
    /* скорость контактного разрыва */
    v_cont = 0.5 * (vl + vr + fr - fl);
}

/* Определение начального приближения для расчета давления на контактном разрыве */
double godunov::Pressure_initial_guess(double pl, double vl, double rhol,
    double cl, double pr, double vr, double rhor, double cr) {

    /* начальное приближение, рассчитанное на освановании рассмотрения линеаризованной системы
       в примитивных переменных */
    double p_lin;
    double p_min, p_max;                /* минимальное и максимальное давления слева и справа от разрыва */
    double p_ratio;                     /* перепад по давлению слева и справа от разрыва */
    double p_ratio_max = 2.0;           /* максимальный перепад по давлению слева и справа от разрыва */
    double p1, p2;              /* вспомогательные переменные для промежуточных расчетов */
    double eps = 1.e-8;                 /* малый эпсилон для критерия сходимости итераций и сравнения вещественных чисел */

    /* Начальное приближение из линейной задачи */
    p_lin = std::max(0.0, 0.5 * (pl + pr) - 0.125 * (vr - vl) * (rhol + rhor) * (cl + cr));
    p_min = std::min(pl, pr);
    p_max = std::max(pl, pr);
    p_ratio = p_max / p_min;

    double g1 = 0.5 * (gamma - 1.0) / gamma;
    double g2 = 0.5 * (gamma + 1.0) / gamma;
    double g3 = 2.0 * gamma / (gamma - 1.0);
    double g4 = 2.0 / (gamma - 1.0);
    double g5 = 2.0 / (gamma + 1.0);
    double g6 = (gamma - 1.0) / (gamma + 1.0);
    double g7 = 0.5 * (gamma - 1.0);


    if ((p_ratio <= p_ratio_max) &&
        ((p_min < p_lin && p_lin < p_max) || (fabs(p_min - p_lin) < eps || fabs(p_max - p_lin) < eps))) {
        /* Начальное приближение из линеаризованной задачи */
        return p_lin;
    }
    else {
        if (p_lin < p_min) {
            /* Начальное приближение по двум волнам разрежения  */

            double PQ = pow(pl / pr, g1);
            double vc = (PQ * vl / cl + vr / cr + g4 * (PQ - 1.0f) / (PQ / cl + 1. / cr));
            double PTL = 1.0 + g7 * (vl - vc) / cl;
            double PTR = 1.0 + g7 * (vc - vr / cr);
            return 0.5 * (pl * pow(PTL, g3) + pr * pow(PTR, g3));

             //pow(((cl + cr - 0.5 * (gamma - 1.0) * (vr - vl)) / (cl / pow(pl, g1) + cr / pow(pr, g1))), 1.0 / g1);
        }
        else {
            /* Начальное приближение по двум ударным волнам  */

            double GEL = sqrt((g5 / rhol) / (g6 * pl + p_lin));
            double GER = sqrt((g5 / rhor) / (g6 * pr + p_lin));
            return (GEL * pl + GER * pr - (vr - vl)) / (GEL + GER);
        }
    }

}

void godunov::Calc_F_and_DF(double curr_press, double p_F_and_DF, double v_F_and_DF, double rho, double c_F_and_DF, double* F, double* DF) {
    double p_ratio, fg, q;          /* вспомогательные переменные */
    p_ratio = curr_press / p_F_and_DF;
    if (curr_press <= p_F_and_DF) {
        /* волна разрежения */
        fg = 2.0 / (gamma - 1.0);
        *F = fg * c_F_and_DF * (pow(p_ratio, 1.0 / fg / gamma) - 1.0);
        *DF = (1.0 / rho / c_F_and_DF) * pow(p_ratio, -0.5 * (gamma + 1.0) / gamma);
    }
    else {
        /* ударная волна */
        q = sqrt(0.5 * (gamma + 1.0) / gamma * p_ratio + 0.5 * (gamma - 1.0) / gamma);
        *F = (curr_press - p_F_and_DF) / c_F_and_DF / rho / q;
        *DF = 0.25 * ((gamma + 1.0) * p_ratio + 3 * gamma - 1.0) / gamma / rho / c_F_and_DF / pow(q, 3.0);
    }
}

void godunov::BoundaryWallNonCons(bool type, Grid& p, Grid& v, Grid& rho) {
	// 0 - левая 1 - правая
   
    if (!type) {
        for (int i = 1; i <= fict; ++i) {

         /*   (p)[fict - i] = (p)[fict + i];
            (v)[fict - i] = -(v)[fict + i];
            (rho)[fict - i] = (rho)[fict + i];*/

            (p)[fict - i] = (p)[fict + i - 1];
            (v)[fict - i] = -(v)[fict + i - 1];
            (rho)[fict - i] = (rho)[fict + i - 1];
           
        }
    }
    else {
        for (int i = 1; i <= fict; ++i) {

           /* (p)[N_x - 1 - fict + i] = (p)[N_x - fict - i - 1];
            (v)[N_x - 1 - fict + i] = -(v)[N_x - fict - i - 1];
            (rho)[N_x - 1 - fict + i] = (rho)[N_x - fict - i - 1];*/

            (p)[N_x - 1 - fict + i] = (p)[N_x - fict - i];
            (v)[N_x - 1 - fict + i] = -(v)[N_x - fict - i];
            (rho)[N_x - 1 - fict + i] = (rho)[N_x - fict - i];
        }
    }

    //if (!type) { // левая 
	//	p_boundary = (*p)[0 + fict];
	//	v_boundary = -(*v)[0 + fict];
	//	rho_boundary = (*rho)[0 + fict];
	//} else { // правая
	//	p_boundary = (*p)[N_x - fict - 1];
	//	v_boundary = -(*v)[N_x - fict - 1];
	//	rho_boundary = (*rho)[N_x - fict - 1];
	//}
	return;
}

void godunov::BoundaryWallCons(bool type, Grid& m, Grid& imp, Grid& e) {
    // 0 - левая 1 - правая
    
    if (!type) {
        
        
        
        for (int i = 1; i <= fict; ++i) {

            //(imp)[fict - i] = -(imp)[fict + i];
            //(m)[fict - i] = (m)[fict + i];
            //(e)[fict - i] = (e)[fict + i];

            (imp)[fict - i] = -(imp)[fict + i - 1];
            (m)[fict - i] = (m)[fict + i - 1];
            (e)[fict - i] = (e)[fict + i - 1];
            
        }
    }
    else {
        for (int i = 1; i <= fict; ++i) {

            /*(imp)[N_x - 1 - fict + i] = -(imp)[N_x - fict - i - 1];
            (m)[N_x - 1 - fict + i] = (m)[N_x - fict - i - 1];
            (e)[N_x - 1 - fict + i] = (e)[N_x - fict - i - 1];*/

            (imp)[N_x - 1 - fict + i] = -(imp)[N_x - fict - i];
            (m)[N_x - 1 - fict + i] = (m)[N_x - fict - i];
            (e)[N_x - 1 - fict + i] = (e)[N_x - fict - i];

        }
    }
    
    //if (!type) { // левая 
    //    m_boundary = (*m)[0 + fict];
    //    imp_boundary = -(*imp)[0 + fict];
    //    e_boundary = (*e)[0 + fict];
    //}
    //else { // правая
    //    m_boundary = (*m)[N_x - fict - 1];
    //    imp_boundary = -(*imp)[N_x - fict - 1];
    //    e_boundary = (*e)[N_x - fict - 1];
    //}
    return;
}

void godunov::BoundaryOpenNonCons(bool type, Grid& p, Grid& v, Grid& rho) {
    // 0 - левая 1 - правая
    
    if (!type) {
        for (int i = 1; i <= fict; ++i) {

         /*   (p)[fict - i] = (p)[fict + i];
            (v)[fict - i] = (v)[fict + i];
            (rho)[fict - i] = (rho)[fict + i];*/

            (p)[fict - i] = (p)[fict + i - 1];
            (v)[fict - i] = (v)[fict + i - 1];
            (rho)[fict - i] = (rho)[fict + i - 1];

        }

    }
    else {

        for (int i = 1; i <= fict; ++i) {

          /*  (p)[N_x - 1 - fict + i] = (p)[N_x - fict - i - 1];
            (v)[N_x - 1 - fict + i] = (v)[N_x - fict - i - 1];
            (rho)[N_x - 1 - fict + i] = (rho)[N_x - fict - i - 1];*/

            (p)[N_x - 1 - fict - i] = (p)[N_x - fict - i];
            (v)[N_x - 1 - fict - i] = (v)[N_x - fict - i];
            (rho)[N_x - 1 - fict - i] = (rho)[N_x - fict - i];
           
        }

    }
    
    
    
    //if (!type) { // левая 
    //    p_boundary = (*p)[0 + fict];
    //    v_boundary = (*vb)[0 + fict];
    //    rho_boundary = (*rho)[0 + fict];
    //}
    //else { // правая
    //    p_boundary = (*p)[N_x - fict - 1];
    //    v_boundary = (*vb)[N_x - fict - 1];
    //    rho_boundary = (*rho)[N_x - fict - 1];
    //}
    return;
}

void godunov::BoundaryOpenCons(bool type, Grid& m, Grid& imp, Grid& e) {
    // 0 - левая 1 - правая
    
    if (!type) {
        for (int i = 1; i <= fict; ++i) {
            /*(imp)[fict - i] = (imp)[fict + i];
            (m)[fict - i] = (m)[fict + i];
            (e)[fict - i] = (e)[fict + i];*/

            (imp)[fict - i] = (imp)[fict + i - 1];
            (m)[fict - i] = (m)[fict + i - 1];
            (e)[fict - i] = (e)[fict + i - 1];
            
        }
    }
    else {
        for (int i = 1; i <= fict; ++i) {
            
           /* (imp)[N_x - 1 - fict + i] = (imp)[N_x - fict - i - 1];
            (m)[N_x - 1 - fict + i] = (m)[N_x - fict - i - 1];
            (e)[N_x - 1 - fict + i] = (e)[N_x - fict - i - 1];*/

            (imp)[N_x - 1 - fict + i] = (imp)[N_x - fict - i];
            (m)[N_x - 1 - fict + i] = (m)[N_x - fict - i];
            (e)[N_x - 1 - fict + i] = (e)[N_x - fict - i];
        }
    }
    
    
    //if (!type) { // левая 
    //    m_boundary = (*m)[0 + fict];
    //    imp_boundary = (*imp)[0 + fict];
    //    e_boundary = (*e)[0 + fict];
    //}
    //else { // правая
    //    m_boundary = (*m)[N_x - fict - 1];
    //    imp_boundary = (*imp)[N_x - fict - 1];
    //    e_boundary = (*e)[N_x - fict - 1];
    //}
    return;
}


void godunov::writeResults(string file_path_res, string file_path_an, string format = ".csv", char _delim = '/t', double _time = 0., int test = 1) {
    ofstream file;
    int cutting_format = 3;
    string s = to_string(ceil(pow(10, cutting_format) * _time) / pow(10, cutting_format));
    s.erase(s.begin() + 2 + cutting_format, s.end());

    file.open(file_path_res + s + format);
    if (file.is_open()) {
        file << _time << endl;

        file << "x coordinate" << _delim << "x component of velocity" << _delim << "density" << _delim << "inner energy" << _delim << "pressure" << endl;

        for (size_t i = fict; i < N_x - fict; i++) {
            file << (*x)[i] << _delim << (*v)[i] << _delim << (*rho)[i] << _delim << (*e)[i] << _delim << (*p)[i] << endl;
        }
    }
    else {
        cerr << "matvey durak" << endl;
        exit(1);
    }
    file.close();

    ofstream file2;
    if (test == 0 || test == 1) {
        file2.open(file_path_an + s + format);
        if (file2.is_open()) {
            file2 << _time << endl;

            file2 << "x coordinate" << _delim << "x component of velocity" << _delim << "density" << _delim << "inner energy" << _delim << "pressure" << endl;

            for (size_t i = fict; i < N_x - fict; i++) {
                file2 << (*x)[i] << _delim << (*v_an)[i] << _delim << (*rho_an)[i] << _delim << (*e_an)[i] << _delim << (*p_an)[i] << endl;
            }
        }
        else {
            cerr << "matvey durak" << endl;
            exit(1);
        }
        file2.close();
    }
    
}

/* l - координата разрыва в задаче Сода */
void godunov::Analytic(double time, double l, double x, double m_l, double impl, double e_l, double m_r, double impr, double e_r, double& p_out, 
    double& v_out, double& rho_out, double& m_out, double& imp_out, double& e_out) {
    
    double cr, cl;
    double p_cont, v_cont, r_cont;
    double p1, v1, r1, p2, v2, r2, p3;
    double e, D, XRW1, XRW2, X_cont, XSW;
    double pl, vl, rl, pr, vr, rr;

    Convert_cons_to_noncons(pl, vl, rl, m_l, impl, e_l);
    Convert_cons_to_noncons(pr, vr, rr, m_r, impr, e_r);
    cl = Sound_velocity(pl, rl);
    cr = Sound_velocity(pr, rr);

    if (2.0 * (cl + cr) / (gamma - 1.0) <= vr - vl) {
        /* случай возникновения вакуума */
        printf("\nContact_pressure_velocity -> vacuum is generated\n");
    }
    /* итерационная процедура расчета давления и скорости газа на контактном разрыве*/
    Contact_pressure_velocity(pl, vl, rl, cl, pr, vr, rr, cr, p_cont, v_cont);
    r_cont = rr * ((gamma - 1) * pr + (gamma + 1) * p_cont) / ((gamma + 1) * pr + (gamma - 1) * p_cont);
    r2 = rl * pow(p_cont / pl, 1 / gamma);
    v1 = (2 / (gamma + 1)) * (cl - (l - x) / time);
    r1 = rl * pow((1 - (gamma - 1) / 2 * v1 / cl), 2 / (gamma - 1));
    p1 = pl * pow((1 - (gamma - 1) / 2 * v1 / cl), 2 * gamma / (gamma - 1));
    D = (p_cont - pr) / (rr * v_cont);
    XRW1 = l - cl * time;
    XRW2 = l - (cl - (gamma + 1) * v_cont / 2) * time;
    X_cont = l + v_cont * time;

    XSW = l + D * time;
    if (x < XRW1) {
        p_out = pl;
        v_out = vl;
        rho_out = rl;
    }
    else if (x < XRW2) {
        p_out = p1;
        v_out = v1;
        rho_out = r1;
    }
    else if (x < X_cont) {
        p_out = p_cont;
        v_out = v_cont;
        rho_out = r2;
    }
    else if (x < XSW) {
        p_out = p_cont;
        v_out = v_cont;
        rho_out = r_cont;
    }
    else {
        p_out = pr;
        v_out = vr;
        rho_out = rr;
    }

    Convert_noncons_to_cons(p_out, v_out, rho_out, m_out, imp_out, e_out);
}

double godunov::min_mod_1972 (double a, double b) {
	if (abs(a) < abs (b)) {
		return a;
	} else {
		return b;
	}
}

double godunov::min_mod_1975 (double a, double b) {
	double c = (a+b)/2.;
	if ((abs(a) <= abs (b)) && (abs(a) <= abs (c))) {
		return a;
	} if ((abs(b) <= abs (a)) && (abs(b) <= abs (c))) {
		return b;
	} else {
		return c;
	}
}

double godunov::min_mod_1984 (double a, double b) {
	if (a*a < a*b) {
		return a;
	} else if (b*b <= a*b) {
		return b;
	} else {
		return 0;
	}
}

void godunov::calc_ei() {
    for(int i = 0; i < N_x; ++i){
        (*ei)[i] = (*p)[i] / ((gamma - 1.) * (*rho)[i]);
    }
}

void godunov::Solver_godunov_kolgan(int test, string file_path_res, string file_path_an, string format, char _delim, int fo, int mod_type) {
	double imp_r = 0.;  // импульс справа
	double imp_l = 0.; // импульс слева
	double p_r = 0.; // давление справа
	double p_l = 0.; // давление слева
	double v_r = 0.; // скорость справа
	double v_l = 0.; // скорость слева
	double m_l = 0., e_l = 0., m_r = 0., e_r = 0.; // масса и полныая энергия
	double rho_l = 0., rho_r = 0.; //плотности слева и справа
	double dt, dx;
	
	double F_m_l = 0., F_imp_l = 0., F_e_l = 0., F_m_r = 0., F_imp_r = 0., F_e_r = 0.; // потоки консервативные
	
	double time = 0.0;
	
	int step_number = 0;
	
	makeGrid();

    Grid* p_an_start = new Grid(N_x);
    Grid* v_an_start = new Grid(N_x);
    Grid* rho_an_start = new Grid(N_x);
    Grid* m_an_start = new Grid(N_x);
    Grid* imp_an_start = new Grid(N_x);
    Grid* e_an_start = new Grid(N_x);
	
    if (mod_type == 1 || mod_type == 2 || mod_type == 3) {
        fict = 2;
    }
    else {
        fict = 1;
    }

	switch(test) {
        case 0:
            Test0(N_x, p, rho, v);
            Test0(N_x, p_an_start, rho_an_start, v_an_start);
            break;
        case 1:
            Test1(N_x, p, rho, v);
            Test1(N_x, p_an_start, rho_an_start, v_an_start);
            break;
        case 2:
            Test2(N_x, p, rho, v);
            Test2(N_x, p_an_start, rho_an_start, v_an_start);
            break;
        case 3:
            Test3(N_x, p, rho, v);
            Test3(N_x, p_an_start, rho_an_start, v_an_start);
            break;
        case 4:
            Test4(N_x, p, rho, v);
            Test4(N_x, p_an_start, rho_an_start, v_an_start);
            break;
        case 5:
            Test5(N_x, p, rho, v);
            Test5(N_x, p_an_start, rho_an_start, v_an_start);
            break;
        default:
            std::cout << "Invalid test number" << std::endl;
            exit(1);
            break;
    }
	
    Convert_noncons_to_cons(p, v, rho, m, imp, e);
    Convert_noncons_to_cons(p_an_start, v_an_start, rho_an_start, m_an_start, imp_an_start, e_an_start);
	
	writeResults(file_path_res, file_path_an, format, _delim, time, test);
	
	double end_time = T_1-T_0;
	
	while (time < end_time) {
		dt = Get_dt();
		for (int i = fict; i < N_x-fict; ++i) {
			Kolgan_left((*m)[i], (*m)[i-1], (*m)[i-2], (*m)[i+1], (*x)[i], (*x)[i-1], (*x)[i-2], (*x)[i+1], m_l, m_r, mod_type);
			Kolgan_left((*imp)[i], (*imp)[i-1], (*imp)[i-2], (*imp)[i+1], (*x)[i], (*x)[i-1], (*x)[i-2], (*x)[i+1], imp_l, imp_r, mod_type);
			Kolgan_left((*e)[i], (*e)[i-1], (*e)[i-2], (*e)[i+1], (*x)[i], (*x)[i-1], (*x)[i-2], (*x)[i+1], e_l, e_r, mod_type);

			
			Kolgan_right((*m)[i], (*m)[i-1], (*m)[i+1], (*m)[i+2], (*x)[i], (*x)[i-1], (*x)[i+1], (*x)[i+2], m_l, m_r, mod_type);
			Kolgan_right((*imp)[i], (*imp)[i-1], (*imp)[i+1], (*imp)[i+2], (*x)[i], (*x)[i-1], (*x)[i+1], (*x)[i+2], imp_l, imp_r, mod_type);
			Kolgan_right((*e)[i], (*e)[i-1], (*e)[i+1], (*e)[i+2], (*x)[i], (*x)[i-1], (*x)[i+1], (*x)[i+2], e_l, e_r, mod_type);
			
			//cout << "pravie: i = " << i << " m_l = " << m_l << " m_r = " << m_r << endl;
			
			Godunov_flux(m_l, imp_l, e_l, m_r, imp_r, e_r, F_m_r, F_imp_r, F_e_r);
			
			
			dx = (Lx_1 - Lx_0)/(N_x);
			
			(*m_new)[i] = (*m)[i] - dt * (F_m_r - F_m_l)/dx;
			(*imp_new)[i] = (*imp)[i] - dt * (F_imp_r - F_imp_l)/dx;
			(*e_new)[i] = (*e)[i] - dt * (F_e_r - F_e_l)/dx;
		}
			
		for (int i = 0; i < N_x; ++i) {
			(*imp)[i] = (*imp_new)[i];
			(*m)[i] = (*m_new)[i];
			(*e)[i] = (*e_new)[i];
			Convert_cons_to_noncons ((*p)[i], (*v)[i], (*rho)[i], (*m)[i], (*imp)[i], (*e)[i]);
		}
		
		for (int i = 0; i <= fict; ++i) {
			(*imp)[i] = (*imp)[fict];
			(*m)[i] = (*m)[fict];
			(*e)[i] = (*e)[fict];
			(*imp)[N_x - 1 - i] = (*imp)[N_x -fict - 1];
			(*m)[N_x - 1 - i] = (*m)[N_x - fict - 1];
			(*e)[N_x - 1 - i] = (*e)[N_x - fict - 1];
		}
		
		for (int i = 0; i < N_x; i++) {
			if (isnan((*m)[i]) || isinf ((*m)[i]))
			cout << "step_number = " << step_number << " i = " << i << " m[" << i << "] = " << (*m)[i] << endl;
		}
		
		time += dt;
		step_number++;
        if (step_number % fo == 0) {
            double l_coord = (*x)[N_x / 2];
			
            for (size_t i = 0; i < N_x; i++) {
                Analytic(time, l_coord, (*x)[i], (*m_an_start)[0], (*imp_an_start)[0], (*e_an_start)[0], (*m_an_start)[N_x - 1], 
                     (*imp_an_start)[N_x - 1], (*e_an_start)[N_x - 1], (*p_an)[i], (*v_an)[i], (*rho_an)[i], (*m_an)[i], (*imp_an)[i], (*e_an)[i]);
                
            }	
            writeResults(file_path_res, file_path_an, format, _delim, time, test);
        }
	}
    delete p_an_start;
    delete v_an_start;
    delete rho_an_start;
    delete m_an_start;
    delete imp_an_start;
    delete e_an_start;
}

//Реконструкция слева
void godunov::Kolgan_left(double qi, double qi_1, double qi_2, double qi1, double xi, double xi_1, double xi_2, double xi1, double& ql, double& qr, int mod_type){
    double a, b, c; //linear coeffs
    
    a = (qi - qi_1)/(xi - xi_1);
    b = (qi1 - qi)/(xi1 - xi);
	
    switch(mod_type){
        case 1:
            c = min_mod_1972(a, b);
            break;
        case 2:
            c = min_mod_1975(a, b);
			break;
        case 3:
            c = min_mod_1984(a, b);
            break;
        default:
            c = min_mod_1984(a, b);
            break;
    }
	if (isinf(c)||isnan(c)) {
		cout << "cl_1 = "<< c << endl;
	}
    qr = -c * (xi - xi_1) / 2. + qi;
    //b = a;
    b = (qi_1 - qi_2)/(xi_1 - xi_2);
    switch(mod_type){
        case 1:
            c = min_mod_1972(a, b);
            break;
        case 2:
            c = min_mod_1975(a, b);
            break;
        case 3:
            c = min_mod_1984(a, b);
            break;
        default:
            c = min_mod_1984(a, b);
            break;
    }
	
	if (isinf(c)||isnan(c)) {
		cout << "cl_2 = "<< c << endl;
	}
    ql = c * (xi - xi_1) / 2. + qi_1;    
}

//Реконструкция справа
void godunov::Kolgan_right(double qi, double qi_1, double qi1, double qi2, double xi, double xi_1, double xi1, double xi2, double& ql, double& qr, int mod_type) {
    double a, b, c; //linear coeffs
    
    a = (qi - qi_1)/(xi - xi_1);
    b = (qi1 - qi)/(xi1 - xi);
	
    switch(mod_type) {
        case 1:
            c = min_mod_1972(a, b);
            break;
        case 2:
            c = min_mod_1975(a, b);
            break;
        case 3:
            c = min_mod_1984(a, b);
            break;
        default:
            c = min_mod_1984(a, b);
            break;
    }
	
    ql = c * (xi1 - xi) / 2. + qi;
	//b=a;
    a = (qi2 - qi1)/(xi2 - xi1);
	
    switch(mod_type) {
        case 1:
            c = min_mod_1972(a, b);
            break;
        case 2:
            c = min_mod_1975(a, b);
            break;
        case 3:
            c = min_mod_1984(a, b);
            break;
        default:
            c = min_mod_1984(a, b);
            break;
    }
    qr = -c * (xi1 - xi) / 2. + qi1;    
}
//                    .-=====-.
//                    | .""". |
//                    | |   | |
//                    | |   | |
//                    | '---' |
//                    |       |
//                    |       |
// .-================-'       '-================-.
//j|  _                                          |
//g|.'o\                                    __   |
//s| '-.'.                                .'o.`  |
// '-==='.'.=========-.       .-========.'.-'===-'
//        '.`'._    .===,     |     _.-' /
//          '._ '-./  ,//\   _| _.-'  _.'
//             '-.| ,//'  \-'  `   .-'
//                `//'_`--;    ;.-'
//                  `\._ ;|    |
//                     \`-'  . |
//                     |_.-'-._|
//                     \  _'_  /
//                     |; -:- | 
//                     || -.- \ 
//                     |;     .\
//                     / `'\'`\-;
//                    ;`   '. `_/
//                    |,`-._;  .;
//                    `;\  `.-'-;
//                     | \   \  |
//                     |  `\  \ |
//                     |   )  | |
//                     |  /  /` /
//                     | |  /|  |
//                     | | / | /
//                     | / |/ /|
//                     | \ / / |
//                     |  /o | |
//                     |  |_/  |
//                     |       |
//                     |       |
//                     |       |
//                     |       |
//                     |       |
//                     |       |
//                     |       |
//                     '-=====-'

void godunov::Solver_godunov_WENO(int test, string file_path_res, string file_path_an, string format, char _delim, int fo) {
   
	double dt, dx;
    //для расщепления потоков WENO
    Grid Fm_plus(N_x), Fimp_plus(N_x), Fe_plus(N_x);
    Grid Fm_minus(N_x), Fimp_minus(N_x), Fe_minus(N_x);
    Grid Fm(N_x, 0.);
    Grid Fimp(N_x, 0.);
    Grid Fe(N_x, 0.);

    //для смещения в РК4
    Grid Flowm_l(N_x), Flowimp_l(N_x), Flowe_l(N_x);
    Grid Flowm_r(N_x), Flowimp_r(N_x), Flowe_r(N_x);
    Grid L_m(N_x), L_imp(N_x), L_e(N_x);

	double time = 0.0;
	
	int step_number = 0;
	
	makeGrid();

    Grid* p_an_start = new Grid(N_x);
    Grid* v_an_start = new Grid(N_x);
    Grid* rho_an_start = new Grid(N_x);
    Grid* m_an_start = new Grid(N_x);
    Grid* imp_an_start = new Grid(N_x);
    Grid* e_an_start = new Grid(N_x);
	

	switch(test) {
        case 0:
            Test0(N_x, p, rho, v);
           /* Test0(N_x, p_left, rho_left, v_left);
            Test0(N_x, p_right, rho_right, v_right);*/
            Test0(N_x, p_an_start, rho_an_start, v_an_start);
            break;
        case 1:
            Test1(N_x, p, rho, v);
            /*Test1(N_x, p_left, rho_left, v_left);
            Test1(N_x, p_right, rho_right, v_right);*/
            Test1(N_x, p_an_start, rho_an_start, v_an_start);
            break;
        case 2:
            Test2(N_x, p, rho, v);
            /*Test2(N_x, p_left, rho_left, v_left);
            Test2(N_x, p_right, rho_right, v_right);*/
       
            break;
        case 3:
            Test3(N_x, p, rho, v);
            /*Test3(N_x, p_left, rho_left, v_left);
            Test3(N_x, p_right, rho_right, v_right);*/
           
            break;
        case 4:
            Test4(N_x, p, rho, v);
            /*Test4(N_x, p_left, rho_left, v_left);
            Test4(N_x, p_right, rho_right, v_right);*/
            
            break;
        case 5:
            Test5(N_x, p, rho, v);
            /*Test5(N_x, p_left, rho_left, v_left);
            Test5(N_x, p_right, rho_right, v_right);*/
            break;
        case 6:
            Test6(N_x, p, rho, v);
           /* Test6(N_x, p_left, rho_left, v_left);
            Test6(N_x, p_right, rho_right, v_right);*/
            break;

        default:
            std::cout << "Invalid test number" << std::endl;
            exit(1);
            break;
    }
	
    Convert_noncons_to_cons(p, v, rho, m, imp, e);
    Convert_noncons_to_cons(p_an_start, v_an_start, rho_an_start, m_an_start, imp_an_start, e_an_start);
	
	writeResults(file_path_res, file_path_an, format, _delim, time, test);
	
	double end_time = T_1-T_0;

    vector<vector<double>> rk = { {1.0,    0.0,       0.0,       0.0},
                                      {1.0,    0.0,       0.0,       1.0},
                                      {0.75, 0.25,       0.0,    0.25},
                                      {1.0 / 3.0,    0.0, 2.0 / 3.0, 2.0 / 3.0} };
    
    

	while (time < end_time) {

        BoundaryWallCons(0, *m, *imp, *e);

        BoundaryOpenCons(1, *m, *imp, *e);

        vector<Grid> m_temp(4, Grid(N_x, 0.));
        vector<Grid> imp_temp(4, Grid(N_x, 0.));
        vector<Grid> e_temp(4, Grid(N_x, 0.));

        m_temp[0] = *m;
        imp_temp[0] = *imp;
        e_temp[0] = *e;

        

        dt = Get_dt();

        for (size_t s = 1; s < 4; ++s) {

            double p_t = 0., v_t = 0., rho_t = 0.;
            
            //считаем стандарртные потоки консервативных переменных
            for (size_t i = 0; i < N_x; ++i) {
                Convert_cons_to_noncons(p_t, v_t, rho_t, m_temp[s - 1][i], imp_temp[s - 1][i], e_temp[s - 1][i]);
                
                Fm[i] = v_t * rho_t;
                Fimp[i] = v_t * v_t * rho_t + p_t;
                Fe[i] = v_t * (e_temp[s - 1][i] + p_t);

                //Diff_flux_ncons_x(p_t, v_t, rho_t, Fm[i], Fimp[i], Fe[i]);
            }

            //считаем максимальную скорость
            double A = 0.;
            for (size_t i = fict; i < N_x - fict; ++i) {
                Convert_cons_to_noncons(p_t, v_t, rho_t, m_temp[s - 1][i], imp_temp[s - 1][i], e_temp[s - 1][i]);
                double c = Sound_velocity(p_t, rho_t);
                if (A < c) {
                    A = c;
                }
            }

            //расчёт расщеплённых потоков
            for (size_t i = fict; i < N_x - fict; ++i) {
                Convert_cons_to_noncons(p_t, v_t, rho_t, m_temp[s - 1][i], imp_temp[s - 1][i], e_temp[s - 1][i]);
                WENO_flux(p_t, v_t, rho_t, A, (Fm)[i], (Fimp)[i], (Fe)[i], (Fm_minus)[i], (Fimp_minus)[i], (Fe_minus)[i], (Fm_plus)[i], (Fimp_plus)[i], (Fe_plus)[i]);
            }

            //реконструкция решения методом WENO
            //забыли про переменные, реконструируем потоки

            WENO_reconstructLeft(Fm_minus, Flowm_l);
            WENO_reconstructRight(Fm_plus, Flowm_r);

            WENO_reconstructLeft(Fimp_minus, Flowimp_l);
            WENO_reconstructRight(Fimp_plus, Flowimp_r);

            WENO_reconstructLeft(Fe_minus, Flowe_l);
            WENO_reconstructRight(Fe_plus, Flowe_r);

            dx = (Lx_1 - Lx_0) / (N_x);

            BoundaryWallCons(0, Flowm_l, Flowimp_l, Flowe_l);
            BoundaryWallCons(0, Flowm_r, Flowimp_r, Flowe_r);
            
            BoundaryOpenCons(1, Flowm_l, Flowimp_l, Flowe_l);
            BoundaryOpenCons(1, Flowm_r, Flowimp_r, Flowe_r);

            //считаем смещение для РК
            for (size_t i = fict + 1; i < N_x - fict - 1; ++i) {
                L_m[i] = ((Flowm_r[i] + Flowm_l[i + 1]) - (Flowm_r[i - 1] + Flowm_l[i])) / dx;

                L_imp[i] = ((Flowimp_r[i] + Flowimp_l[i + 1]) - (Flowimp_r[i - 1] + Flowimp_l[i])) / dx;

                L_e[i] = ((Flowe_r[i] + Flowe_l[i + 1]) - (Flowe_r[i - 1] + Flowe_l[i])) / dx;
            }

            for (size_t i = fict; i < N_x - fict; ++i) {

                m_temp[s][i] = rk[s][0] * m_temp[0][i]
                    + rk[s][1] * m_temp[1][i]
                    + rk[s][2] * m_temp[2][i]
                    - rk[s][3] * dt * L_m[i];

                imp_temp[s][i] = rk[s][0] * imp_temp[0][i]
                    + rk[s][1] * imp_temp[1][i]
                    + rk[s][2] * imp_temp[2][i]
                    - rk[s][3] * dt * L_imp[i];

                e_temp[s][i] = rk[s][0] * e_temp[0][i]
                    + rk[s][1] * e_temp[1][i]
                    + rk[s][2] * e_temp[2][i]
                    - rk[s][3] * dt * L_e[i];
            }
           
           

        }
        
        for (size_t i = fict; i < N_x - fict; ++i) {

            (*m)[i] =  m_temp[3][i];

            (*imp)[i] = imp_temp[3][i];

            (*e)[i] = e_temp[3][i];

            if ((*m)[i] < 0.) {
                cout << "t = " << time << ", i = " << i << ", m = " << (*m)[i] << endl;
            }
            else if ((*e)[i] < 0.) {
                cout << "t = " << time << ", i = " << i << ", e = " << (*e)[i] << endl;
            }
        }

        BoundaryWallCons(0, *m, *imp, *e);

        BoundaryOpenCons(1, *m, *imp, *e);

        Convert_cons_to_noncons(p, v, rho, m, imp,e);
		
		time += dt;
		step_number++;
        if (step_number % fo == 0) {
            double l_coord = (*x)[N_x / 2];
            if (test == 0 || test == 1) {
                for (size_t i = 0; i < N_x; i++) {
                                Analytic(time, l_coord, (*x)[i], (*m_an_start)[0], (*imp_an_start)[0], (*e_an_start)[0], (*m_an_start)[N_x - 1], 
                                     (*imp_an_start)[N_x - 1], (*e_an_start)[N_x - 1], (*p_an)[i], (*v_an)[i], (*rho_an)[i], (*m_an)[i], (*imp_an)[i], (*e_an)[i]);
                }
            }
            	
            writeResults(file_path_res, file_path_an, format, _delim, time, test);
        }
	}
    delete p_an_start;
    delete v_an_start;
    delete rho_an_start;
    delete m_an_start;
    delete imp_an_start;
    delete e_an_start;
}

void godunov::WENO_flux(double p, double v, double rho, double c_max, double Fm, double Fimp, double Fe,
    double& Fm_minus, double& Fimp_minus, double& Fe_minus, double& Fm_plus, double& Fimp_plus, double& Fe_plus) {
    
    double m = 0., imp = 0., e = 0.;

    Convert_noncons_to_cons(p, v, rho, m, imp, e);

    //расщепление потоков
    Fm_minus = 0.5 * (Fm - c_max * m);
    Fm_plus = 0.5 * (Fm + c_max * m);

    Fimp_minus = 0.5 * (Fimp - c_max * imp);
    Fimp_plus = 0.5 * (Fimp + c_max * imp);

    Fe_minus = 0.5 * (Fe - c_max * e);
    Fe_plus = 0.5 * (Fe + c_max * e);

    return;
}
void godunov::WENO_reconstructLeft(Grid v, Grid& Flow) {

    double eps = 1e-8;
    vector<double> gamma = { 0.3, 0.6, 0.1 };
    vector<double> a(3);
    vector<double> h(3);
    vector<double> beta(3);
    vector<double> omega(3);

    //потоки для разности L
   

    //задействовано 5 точек
    for (size_t i = fict; i < N_x - fict; ++i) {
         
        //шаблоны: для левого
        
        h[0] = 11.0 / 6.0 * v[i] - 7.0 / 6.0 * v[i + 1] + 1.0 / 3.0 * v[i + 2];

        h[1] = 1.0 / 3.0 * v[i - 1] + 5.0 / 6.0 * v[i] - 1.0 / 6.0 * v[i + 1];

        h[2] = -1.0 / 6.0 * v[i - 2] + 5.0 / 6.0 * v[i - 1] + 1.0 / 3.0 * v[i];

        //индикаторы гладкости
        
        beta[0] = 13.0 / 12.0 * pow(v[i] - 2.0 * v[i + 1] + v[i + 2], 2.0) + 1.0 / 4.0 * pow(3.0 * v[i] - 4.0 * v[i + 1] + v[i + 2], 2.0);

        beta[1] = 13.0 / 12.0 * pow(v[i - 1] - 2.0 * v[i] + v[i + 1], 2.0) + 1.0 / 4.0 * pow(v[i - 1] - v[i + 1], 2.0);

        beta[2] = 13.0 / 12.0 * pow(v[i - 2] - 2.0 * v[i - 1] + v[i], 2.0) + 1.0 / 4.0 * pow(v[i - 2] - 4.0 * v[i - 1] + 3.0 * v[i], 2.0);

        double summ = 0.;

        for (int j = 0; j < 3; ++j) {
            if (beta[j] != 0.) a[j] = gamma[j] / pow(beta[j], 2);

            else  a[j] = gamma[j] / pow(eps, 2.0);
            
            summ += a[j];
        }

        for (int j = 0; j < 3; ++j) omega[j] = a[j] / summ;
        
        //считаем левый "поток" 
       
        Flow[i] = inner_product(omega.begin(), omega.end(), h.begin(), 0.);

    }
    return;
}
void godunov::WENO_reconstructRight(Grid v, Grid& Flow) {

    double eps = 1e-8;
    vector<double> gamma = { 0.3, 0.6, 0.1 };
    vector<double> a(3);
    vector<double> h(3);
    vector<double> beta(3);
    vector<double> omega(3);
    
    //разворачиваем gamma
    reverse(gamma.begin(), gamma.end());

    //потоки для разности L

    //задействовано 5 точек
    for (size_t i = fict; i < N_x - fict; ++i) {

        //шаблоны: для правого
        
        h[0] = 1.0 / 3.0 * v[i] + 5.0 / 6.0 * v[i + 1] - 1.0 / 6.0 * v[i + 2];

        h[1] = -1.0 / 6.0 * v[i - 1] + 5.0 / 6.0 * v[i] + 1.0 / 3.0 * v[i + 1];

        h[2] = 1.0 / 3.0 * v[i - 2] - 7.0 / 6.0 * v[i - 1] + 11.0 / 6.0 * v[i];

        //индикаторы гладкости
        beta[0] = 13.0 / 12.0 * pow(v[i] - 2.0 * v[i + 1] + v[i + 2], 2.0) + 1.0 / 4.0 * pow(3.0 * v[i] - 4.0 * v[i + 1] + v[i + 2], 2.0);

        beta[1] = 13.0 / 12.0 * pow(v[i - 1] - 2.0 * v[i] + v[i + 1], 2.0) + 1.0 / 4.0 * pow(v[i - 1] - v[i + 1], 2.0);

        beta[2] = 13.0 / 12.0 * pow(v[i - 2] - 2.0 * v[i - 1] + v[i], 2.0) + 1.0 / 4.0 * pow(v[i - 2] - 4.0 * v[i - 1] + 3.0 * v[i], 2.0);

        double summ = 0.;

        for (int j = 0; j < 3; ++j) {
            if (beta[j] != 0.) a[j] = gamma[j] / pow(beta[j], 2);

            else  a[j] = gamma[j] / pow(eps, 2.0);

            summ += a[j];
        }

        for (int j = 0; j < 3; ++j) omega[j] = a[j] / summ;

        //считаем правый "поток" 

         Flow[i] = inner_product(omega.begin(), omega.end(), h.begin(), 0.);
    }
    return;
}
void godunov::Solver_mccormak(int test, string file_path_res, string file_path_an, string format, char _delim, int fo) {
    double dt, dx;
    
        double time = 0.0;
   
        int step_number = 0;
    
        vector <double> a = { 0, 1 };
        vector<double> b = { 0.5, 0.5 };

        makeGrid();
        
        Grid* p_an_start = new Grid(N_x);
        Grid* v_an_start = new Grid(N_x);
        Grid* rho_an_start = new Grid(N_x);
        Grid* m_an_start = new Grid(N_x);
        Grid* imp_an_start = new Grid(N_x);
        Grid* e_an_start = new Grid(N_x);
        

        
        switch (test) {
        case 0:
            Test0(N_x, p, rho, v);
            Test0(N_x, p_an_start, rho_an_start, v_an_start);
            break;
        case 1:
            Test1(N_x, p, rho, v);
            Test1(N_x, p_an_start, rho_an_start, v_an_start);
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
        case 6:
            Test6(N_x, p, rho, v);
            break;
        default:
            std::cout << "Invalid test number" << std::endl;
            exit(1);
            break;
        }
        
        Convert_noncons_to_cons(p, v, rho, m, imp, e);
        Convert_noncons_to_cons(p_an_start, v_an_start, rho_an_start, m_an_start, imp_an_start, e_an_start);
        
        writeResults(file_path_res, file_path_an, format, _delim, time, test);
        
        double end_time = T_1 - T_0;
        size_t iter = 0;
        vector<Grid> u = { (*m), (*imp), (*e) };

        while (time < end_time and iter < max_iter) {
            dt = Get_dt();
        
            // потоки консервативные Маккормака
            Grid F_m(N_x + 1, 0.), F_imp(N_x + 1, 0.), F_e(N_x + 1, 0.);
            
            BoundaryOpenCons(0, u[0], u[1], u[2]);
        
            BoundaryWallCons(1, u[0], u[1], u[2]);
        
            for (size_t i = 0; i < N_x; ++i) {
                Convert_cons_to_noncons((*p)[i], (*v)[i], (*rho)[i], u[0][i], u[1][i], u[2][i]);
            }
        
            dx = (Lx_1 - Lx_0) / (N_x);
        
            vector<Grid> Flux = { F_m, F_imp, F_e };

            vector<Grid> u_mon = vector<Grid>(3, Grid(N_x, 0.));
            vector<Grid> u_copy = u;
            vector<Grid> u_after_d;

            //predictor

            //0 - left, 1 - right
            bool leftOrRight = false;
            
        
                for (int i = fict; i < N_x - fict; ++i) {
                    Mccormak_flux(u, Flux, i, leftOrRight);
                    u[0][i] -= Flux[0][i] * dt / dx;
                    u[1][i] -= Flux[1][i] * dt / dx;
                    u[2][i] -= Flux[2][i] * dt / dx;
                }
        
            //corrector
            
            leftOrRight = true;
            

            for (size_t i = fict; i < N_x - fict; ++i) {
                Mccormak_flux(u, Flux, i, leftOrRight);
                u_copy[0][i] = 0.5 * ((*m)[i] + u[0][i]) - 0.5 * Flux[0][i] * dt / dx;
                u_copy[1][i] = 0.5 * ((*imp)[i] + u[1][i]) - 0.5 * Flux[1][i] * dt / dx;
                u_copy[2][i] = 0.5 * ((*e)[i] + u[2][i]) - 0.5 * Flux[2][i] * dt / dx;
                
                (*m)[i] = u_copy[0][i];
                (*imp)[i] = u_copy[1][i];
                (*e)[i] = u_copy[2][i];
            }
        
            u = u_copy;
           
            /*Diffusion(u_copy, u_mon);
        
            for (int i = 0; i < N_x; ++i) {
        
                for (size_t k = 0; k < 3; ++k) {
                    u_mon[k][i] += u[k][i];
                }
            }*/
        
            //technical step to do (1 + N) * u_d
             
            //u_after_d = u + D(u), needed for (1 + N)(1 + D) (u) 
                
           /* u_after_d = u_mon;
            Antidiffusion(u_after_d);
        
            for (int i = fict; i < N_x - fict; ++i) {
        
                for (size_t k = 0; k < 3; ++k) {
                    u_mon[k][i] += u_after_d[k][i];
                }
            }*/
            //monotone sol
            /*u = u_mon;*/
        
        
            time += dt;
            step_number++;
            if (step_number % fo == 0) {
                    
                double l_coord;
        
                if (test == 0) {
                    l_coord = (*x)[fict + 1];
                }
                else  l_coord = (*x)[N_x / 2];
        
                if (test == 0 || test == 1) {
                    for (size_t i = 0; i < N_x; i++) {
                                    Analytic(time, l_coord, (*x)[i], (*m_an_start)[0], (*imp_an_start)[0], (*e_an_start)[0], (*m_an_start)[N_x - 1],
                                        (*imp_an_start)[N_x - 1], (*e_an_start)[N_x - 1], (*p_an)[i], (*v_an)[i], (*rho_an)[i], (*m_an)[i], (*imp_an)[i], (*e_an)[i]);
        
                                }
                }
                    
                writeResults(file_path_res, file_path_an, format, _delim, time, test);
            }
        }
        delete p_an_start;
        delete v_an_start;
        delete rho_an_start;
        delete m_an_start;
        delete imp_an_start;
        delete e_an_start;

}

//hard realization, but iterative. Not sure if it's correct
//void godunov::Solver_mccormak(int test, string file_path_res, string file_path_an, string format, char _delim, int fo) {
//   
//    double dt, dx;
//
//    double time = 0.0;
//
//    int step_number = 0;
//
//    vector <double> a = { 0, 1 };
//    vector<double> b = { 0.5, 0.5 };
//
//    //векторы приращений переменных
//    
//
//    makeGrid();
//
//    Grid* p_an_start = new Grid(N_x);
//    Grid* v_an_start = new Grid(N_x);
//    Grid* rho_an_start = new Grid(N_x);
//    Grid* m_an_start = new Grid(N_x);
//    Grid* imp_an_start = new Grid(N_x);
//    Grid* e_an_start = new Grid(N_x);
//
//
//    switch (test) {
//    case 0:
//        Test0(N_x, p, rho, v);
//        Test0(N_x, p_an_start, rho_an_start, v_an_start);
//        break;
//    case 1:
//        Test1(N_x, p, rho, v);
//        Test1(N_x, p_an_start, rho_an_start, v_an_start);
//        break;
//    case 2:
//        Test2(N_x, p, rho, v);
//        
//        break;
//    case 3:
//        Test3(N_x, p, rho, v);
//        
//        break;
//    case 4:
//        Test4(N_x, p, rho, v);
//        
//        break;
//    case 5:
//        Test5(N_x, p, rho, v);
//        
//        break;
//    case 6:
//        Test6(N_x, p, rho, v);
//        break;
//    default:
//        std::cout << "Invalid test number" << std::endl;
//        exit(1);
//        break;
//    }
//
//    Convert_noncons_to_cons(p, v, rho, m, imp, e);
//    Convert_noncons_to_cons(p_an_start, v_an_start, rho_an_start, m_an_start, imp_an_start, e_an_start);
//
//    writeResults(file_path_res, file_path_an, format, _delim, time, test);
//
//    double end_time = T_1 - T_0;
//    size_t iter = 0;
//    vector<Grid> u = { (*m), (*imp), (*e) };
//    while (time < end_time and iter < max_iter) {
//        dt = Get_dt();
//
//        // потоки консервативные Маккормака
//        Grid F_m(N_x + 1, 0.), F_imp(N_x + 1, 0.), F_e(N_x + 1, 0.);
//        vector<Grid> h(3, Grid(N_x, 0.));
//        vector<Grid> H(3, Grid(N_x, 0.));
//
//        
//        BoundaryOpenCons(0, *m, *imp, *e);
//
//        BoundaryWallCons(1, *m, *imp, *e);
//
//        Convert_cons_to_noncons(p, v, rho, m, imp, e);
//
//        dx = (Lx_1 - Lx_0) / (N_x);
//
//        //Шаг предиктор: считаем потоки по неконс перем
//
//
//        vector<Grid> Flux = { F_m, F_imp, F_e };
//        
//        vector<Grid> u_mon = vector<Grid>(3, Grid(N_x, 0.));
//        vector<Grid> u_copy;
//        vector<Grid> u_after_d;
//        //s = 0 - предиктор, s = 1 - корректор
//        for (size_t s = 0; s < 2; ++s) {
//
//            for (int i = fict + s; i < N_x - fict - s; ++i) {
//
//                for (size_t k = 0; k < 3; ++k) {
//                    u[k][i] += a[s] * h[k][i];
//                }
//            }
//
//            Mccormak_flux(u, Flux);
//
//            for (int i = fict + s; i < N_x - fict + s; ++i) {
//
//                for (size_t k = 0; k < 3; ++k) {
//
//                    if (s % 2 == 0)     h[k][i] = -dt * (Flux[k][i + 1] - Flux[k][i]) / dx;
//                    else            h[k][i] = -dt * (Flux[k][i] - Flux[k][i - 1]) / dx;
//
//                    H[k][i] += b[s] * h[k][i];
//
//                }
//            }
//
//        }
//
//        for (size_t i = fict; i < N_x - fict; ++i) {
//            
//            (*m)[i] += H[0][i];
//            
//            (*imp)[i] += H[1][i];
//
//            (*e)[i] += H[2][i];
//        }
//
//        u_copy = u;
//        Diffusion(u_copy, u_mon);
//
//        for (int i = 0; i < N_x; ++i) {
//
//            for (size_t k = 0; k < 3; ++k) {
//                u_mon[k][i] += u[k][i];
//            }
//        }
//
//        //technical step to do (1 + N) * u_d
//        // 
//        // u_after_d = u + D(u), needed for (1 + N)(1 + D) (u) 
//        
//        u_after_d = u_mon;
//        Antidiffusion(u_copy, u_after_d);
//
//        for (int i = fict; i < N_x - fict; ++i) {
//
//            for (size_t k = 0; k < 3; ++k) {
//                u_mon[k][i] += u_after_d[k][i];
//            }
//        }
//        //monotone sol
//        u = u_mon;
//
//
//        time += dt;
//        step_number++;
//        if (step_number % fo == 0) {
//            
//            double l_coord;
//
//            if (test == 0) {
//                l_coord = (*x)[fict + 1];
//            }
//            else  l_coord = (*x)[N_x / 2];
//
//            if (test == 0 || test == 1) {
//                for (size_t i = 0; i < N_x; i++) {
//                                Analytic(time, l_coord, (*x)[i], (*m_an_start)[0], (*imp_an_start)[0], (*e_an_start)[0], (*m_an_start)[N_x - 1],
//                                    (*imp_an_start)[N_x - 1], (*e_an_start)[N_x - 1], (*p_an)[i], (*v_an)[i], (*rho_an)[i], (*m_an)[i], (*imp_an)[i], (*e_an)[i]);
//
//                            }
//            }
//            
//            writeResults(file_path_res, file_path_an, format, _delim, time, test);
//        }
//    }
//    delete p_an_start;
//    delete v_an_start;
//    delete rho_an_start;
//    delete m_an_start;
//    delete imp_an_start;
//    delete e_an_start;
//}

void godunov::Mccormak_flux(vector<Grid> u, vector<Grid>& Flux, size_t i, bool leftOrRight) {

    double p_pl = 0., v_pl = 0., rho_pl = 0.;
    double p_mn = 0., v_mn = 0., rho_mn = 0.;
    //right
    if (leftOrRight) {

        Convert_cons_to_noncons(p_pl, v_pl, rho_pl, u[0][i + 1], u[1][i + 1], u[2][i + 1]);
        Convert_cons_to_noncons(p_mn, v_mn, rho_mn, u[0][i], u[1][i], u[2][i]);

        Flux[0][i] = rho_pl * v_pl - rho_mn * v_mn;
        Flux[1][i] = rho_pl * v_pl * v_pl - rho_mn * v_mn * v_mn + p_pl - p_mn;
        Flux[2][i] = v_pl * (0.5 * rho_pl * v_pl * v_pl + p_pl / (gamma - 1.) + p_pl) - v_mn * (0.5 * rho_mn * v_mn * v_mn + p_mn / (gamma - 1.) + p_mn);

    }
    //left
    else {
       
        Convert_cons_to_noncons(p_pl, v_pl, rho_pl, u[0][i], u[1][i], u[2][i]);
        Convert_cons_to_noncons(p_mn, v_mn, rho_mn, u[0][i - 1], u[1][i - 1], u[2][i - 1]);
           
        Flux[0][i] = rho_pl * v_pl - rho_mn * v_mn;
        Flux[1][i] = rho_pl * v_pl * v_pl - rho_mn * v_mn * v_mn + p_pl - p_mn;
        Flux[2][i] = v_pl * (0.5 * rho_pl * v_pl * v_pl + p_pl / (gamma - 1.) + p_pl) - v_mn * (0.5 * rho_mn * v_mn * v_mn + p_mn / (gamma - 1.) + p_mn);

    }
    return;
}

void godunov::Diffusion(vector<Grid> u, vector<Grid>& Du) {
    double C_diff = 0.1;

    //Получааем диффузию

    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = fict; i < N_x - fict; ++i) Du[k][i] = C_diff * (u[k][i + 1] - 2 * u[k][i] + u[k][i - 1]);

        Du[k][0] = Du[k][1];
        Du[k][N_x - 1] = Du[k][N_x - 2];
    }
    
}

void godunov::Antidiffusion(vector<Grid>& Nu) {
    double C_diff = 0.1, C_antidiff = 1.0;

    vector<Grid> del_u(3, Grid(N_x));
    vector<Grid> del_u_d(3, Grid(N_x));
    vector<Grid> phi_c(3, Grid(N_x));

    del_u = Nu;

    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = fict; i < N_x - fict; ++i) del_u_d[k][i] = del_u[k][i] - del_u[k][i - 1];

        del_u_d[k][0] = del_u_d[k][1];
    }

    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = fict; i < N_x - fict; ++i) {

            //sign delta_u_i+1/2
            double s = (del_u_d[k][i]) > 0 ? +1. : (del_u_d[k][i]) < 0 ? -1. : 0.;

            //abs(phi_i+1/2) = C_diff * abs(u[k][i + 1] - u[k][i])
            phi_c[k][i] = s * max({ 0., min({ min({ s * del_u_d[k][i - 1], s * del_u_d[k][i + 2]}), abs(C_diff * del_u_d[k][i + 1])})});
        }

    }
    
    //Получааем антидиффузию
    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = fict; i < N_x - fict; ++i) Nu[k][i] = - C_antidiff * (phi_c[k][i] - phi_c[k][i - 1]);

        Nu[k][0] = Nu[k][1];
        Nu[k][N_x - 2] = Nu[k][N_x - 3];
        Nu[k][N_x - 1] = Nu[k][N_x - 2];
    }

}