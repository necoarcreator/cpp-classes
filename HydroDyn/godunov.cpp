#include "godunov.h"

/* ЧТЕНИЕ ВХОДНЫХ ДАННЫХ */

struct Parameters {
    int cells_number_x;           /* число ячеек x*/

    double stop_time;           /* момент времени, для которого строится точное решение */
    std::vector<PVR> params;    /* вектор переменных различных областей */
    std::vector<double> l;      /* вектор разделителей областей */
    double g;                   /* показатель адибаты */
    double CFL;                 /* число Куранта */
    int S_type;                 /* 1 - G, 2 - GK */
    double ax;                   /* начало расчетной области x*/
    double bx;                   /* конец расчетной области x*/

    int fo;                     /* частота вывода файлов */

    double type_b_left;         /* тип левой границы */
    double type_b_right;        /* тип правой границы */


};
void read_parameters(struct Parameters* params) {
    params->g = 1.4;
    /* НАЧАЛЬНЫЕ УСЛОВИЯ */
    params->ax = 0.0;
    params->bx = 1;
    /* ПАРАМЕТРЫ РАСЧЕТА*/
    params->cells_number_x = 100;
    params->stop_time = 2;
    params->CFL = 0.2;
    params->S_type = 1;
    params->fo = 25;
    params->type_b_left = 0;
    params->type_b_right = 0;

}

/* БЛОК РАБОТЫ С ПАМЯТЬЮ */
/* Выделение памяти под одномерный массив элементов */
template <typename T>
void Allocate(int size, T*& mass) {
    mass = new T[size];
}
template <class First, class... Other>
void Allocate(int Nx, First& first, Other&... other) {
    Allocate(Nx, first);
    Allocate(Nx, other...);
}

/* Очищение памяти одномерного массива элементов */
template<class T>
void Delete(T& mass) {
    delete[] mass;
}
template <class First, class... Other>
void Delete(First& first, Other&... other) {
    Delete(first);
    Delete(other...);
}

/* ПЕРЕХОД ОТ КОНСЕРВАТИВНЫХ ПЕРЕМЕННЫХ И ОБРАТНО */
void godunov::Convert_cons_to_noncons(struct Parameters* params, double& p, double& vx, double& r, double& m, double& impx, double& e) {
    double g = params->g;
    p = (g - 1.0) * (e - 0.5 * (pow(impx, 2.0)) / m);
    vx = impx / m;
    r = m;
}
void godunov::Convert_noncons_to_cons(struct Parameters* params, double& p, double& vx, double& r, double& m, double& impx, double& e) {
    double g = params->g;
    m = r;
    impx = r * vx;
    e = 0.5 * r * (pow(vx, 2.0)) + p / (g - 1.0);
}

/* РАСЧЕТ СКОРОСТИ ЗВУКА */
double godunov::Sound_velocity(struct Parameters* params, double p, double r) {
    double g = params->g;
    return std::sqrt(g * p / r);
}

/* ИНИЦИАЛИЗАЦИЯ СЕТКИ */
void godunov::Build_grid(struct Parameters* params, double* xc, double* x) {
    double hx, hy;   /* шаг сетки */
    double right_boundary_x = params->bx;
    double left_boundary_x = params->ax;
    int Nx = params->cells_number_x;
    hx = (right_boundary_x - left_boundary_x) / Nx;
    /* координаты узлов */
    for (int i = 0; i < Nx + 1; ++i) {
        x[i] = left_boundary_x + i * hx;
    }
    /* координаты центров ячеек */
    for (int i = 0; i < Nx; ++i) {
        xc[i] = 0.5 * (x[i] + x[i + 1]);
    }
}

/* ГРАНИЧНЫЕ УСЛОВИЯ */
void godunov::Boundary_x(double p, double vx, double r, double& pb, double& vxb, double& rb, int b_type) {
    /* стенка */
    if (b_type == 0) {
        pb = p;
        vxb = -vx;
        rb = r;
    }
    else  if (b_type == 1) {/* свободная */
        pb = p;
        vxb = vx;
        rb = r;
    }
}


/* Расчет функции F, определяющей скорость газа на контактном разрыве, и ее производной по давлению среды DF
   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 3rd Edition. - Springer,
   1999. - P. 158. - Subroutine PREFUN.
   + поправка на двучленное уравнение состояния: */
void godunov::Calc_F_and_DF(struct Parameters* params, double curr_press, double p, double v, double r, double c, double* F, double* DF) {
    double g = params->g;           /* показатель адиабаты */
    double p_ratio, fg, q;          /* вспомогательные переменные */
    p_ratio = curr_press / p;
    if (curr_press <= p) {
        /* волна разрежения */
        fg = 2.0 / (g - 1.0);
        *F = fg * c * (pow(p_ratio, 1.0 / fg / g) - 1.0);
        *DF = (1.0 / r / c) * pow(p_ratio, -0.5 * (g + 1.0) / g);
    }
    else {
        /* ударная волна */
        q = sqrt(0.5 * (g + 1.0) / g * p_ratio + 0.5 * (g - 1.0) / g);
        *F = (curr_press - p) / c / r / q;
        *DF = 0.25 * ((g + 1.0) * p_ratio + 3 * g - 1.0) / g / r / c / pow(q, 3.0);
    }
}

/* Определение начального приближения для расчета давления на контактном разрыве
   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 157. - Subroutine GUESSP.
   Возвращает искомое начальное приближения */
double godunov::Pressure_initial_guess(struct Parameters* params, double pl, double vl, double rl,
    double cl, double pr, double vr, double rr, double cr) {
    double g = params->g;                           /* показатель адиабаты */

    /* начальное приближение, рассчитанное на освановании рассмотрения линеаризованной системы
       в примитивных переменных */
    double p_lin;
    double p_min, p_max;                /* минимальное и максимальное давления слева и справа от разрыва */
    double p_ratio;                     /* перепад по давлению слева и справа от разрыва */
    double p_ratio_max = 2.0;           /* максимальный перепад по давлению слева и справа от разрыва */
    double p1, p2, g1, g2;              /* вспомогательные переменные для промежуточных расчетов */
    double eps = 1.e-8;                 /* малый эпсилон для критерия сходимости итераций и сравнения вещественных чисел */

    /* Начальное приближение из линейной задачи
       Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
       1999. - P. 128. - Formula (4.47). */
    p_lin = std::max(0.0, 0.5 * (pl + pr) - 0.125 * (vr - vl) * (rl + rr) * (cl + cr));
    p_min = std::min(pl, pr);
    p_max = std::max(pl, pr);
    p_ratio = p_max / p_min;

    if ((p_ratio <= p_ratio_max) &&
        ((p_min < p_lin && p_lin < p_max) || (fabs(p_min - p_lin) < eps || fabs(p_max - p_lin) < eps))) {
        /* Начальное приближение из линеаризованной задачи */
        return p_lin;
    }
    else {
        if (p_lin < p_min) {
            /* Начальное приближение по двум волнам разрежения
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 302. - Formula (9.32) + поправка на двучленное уравнение состояния */
            g1 = 0.5 * (g - 1.0) / g;
            return pow(((cl + cr - 0.5 * (g - 1.0) * (vr - vl)) / (cl / pow(pl, g1) + cr / pow(pr, g1))), 1.0 / g1);
        }
        else {
            /* Начальное приближение по двум ударным волнам
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 128. - Formula (4.48) + поправка на двучленное уравнение состояния */
            g1 = 2.0 / (g + 1.0);
            g2 = (g - 1.0) / (g + 1.0);
            p1 = std::sqrt(g1 / rl / (g2 * pl + p_lin));
            p2 = std::sqrt(g1 / rr / (g2 * pr + p_lin));
            return (p1 * pl + p2 * pr - (vr - vl)) / (p1 + p2);
        }
    }

}

/* Итерационная процедура расчета давления и скорости на контактном разрыве
   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 3rd Edition. - Springer,
   2009. - P. 155. - Subroutine STARPU.*/

void godunov::Contact_pressure_velocity(struct Parameters* params, double pl, double vl, double rl, double cl,
    double pr, double vr, double rr, double cr, double& p_cont, double& v_cont) {
    double p_old;         /* значение давления на предыдущей итерации */
    double fl, fr;        /* значения функций */
    double fld, frd;      /* значения производных */
    int iter_num = 0;     /* количество проведенных итераций */
    int iter_max = 300;   /* максимальное количество итераций */
    double criteria;      /* переменная для определения сходимости */
    double g = params->g; /* показатель адиабаты */
    double eps = 1.e-8;
    if (2.0 * (cl + cr) / (g - 1.0) <= vr - vl) {
        /* случай возникновения вакуума */
        printf("\nContact_pressure_velocity -> vacuum is generated\n");
    }
    /* расчет начального приближения для давления */
    p_old = Pressure_initial_guess(params, pl, vl, rl, cl, pr, vr, rr, cr);
    if (p_old < 0.0) {
        printf("\nContact_pressure_velocity -> initial pressure guess is negative ");
    }
    /* решение нелинейного уравнения для нахождения давления на контактном разрыве методом Ньютона-Рафсона */
    do {
        Calc_F_and_DF(params, p_old, pl, vl, rl, cl, &fl, &fld);
        Calc_F_and_DF(params, p_old, pr, vr, rr, cr, &fr, &frd);
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


/* Функция отбора решения
   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine SAMPLE.
*/
void godunov::Sample_solid_solution(struct Parameters* params, double pl, double vl, double rl, double cl, double pr, double vr, double rr, double cr,
    double p_cont, double v_cont, double s, double& p_res, double& v_res, double& r_res) {

    double g1, g2, g3, g4, g5, g6, g7;      /* вспомогательные переменные, производные от показателя адиабаты,
                                               в соответствии с Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer, 1999. - P. 153. */

                                               /* скорости левых волн */
    double shl, stl;        /* скорости "головы" и "хвоста" левой волны разрежения */
    double sl;              /* скорость левой ударной волны */

    /* скорости правых волн */
    double shr, str;        /* скорости "головы" и "хвоста" правой волны разрежения */
    double sr;              /* скорость правой ударной волны */

    double cml, cmr;        /* скорости звука слева и справа от контактного разрыва */
    double c;               /* локальная скорость звука внутри волны разрежения */
    double p_ratio;
    double r, v, p;         /* отобранные значения объемной доли, плотности, скорости и давления */

    /* производные от показателя адиабаты */
    g1 = 0.5 * (params->g - 1.0) / params->g;
    g2 = 0.5 * (params->g + 1.0) / params->g;
    g3 = 2.0 * params->g / (params->g - 1.0);
    g4 = 2.0 / (params->g - 1.0);
    g5 = 2.0 / (params->g + 1.0);
    g6 = (params->g - 1.0) / (params->g + 1.0);
    g7 = 0.5 * (params->g - 1.0);

    if (s <= v_cont) {
        /* рассматриваемая точка - слева от контактного разрыва */
        if (p_cont <= pl) {
            /* левая волна разрежения */
            shl = vl - cl;
            if (s <= shl) {
                /* параметры слева от разрыва */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                cml = cl * pow(p_cont / pl, g1);
                stl = v_cont - cml;
                if (s > stl) {
                    /* параметры слева от контактного разрыва */
                    r = rl * pow(p_cont / pl, 1.0 / params->g);
                    v = v_cont;
                    p = p_cont;
                }
                else {
                    /* параметры внутри левой волны разрежения */
                    v = g5 * (cl + g7 * vl + s);
                    c = g5 * (cl + g7 * (vl - s));
                    r = rl * pow(c / cl, g4);
                    p = pl * pow(c / cl, g3);
                }
            }
        }
        else {
            /* левая ударная волна */
            p_ratio = p_cont / pl;
            sl = vl - cl * std::sqrt(g2 * p_ratio + g1);
            if (s <= sl) {
                /* параметры слева от разрыва */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                /* параметры за левой ударной волной */
                r = rl * (p_ratio + g6) / (p_ratio * g6 + 1.0);
                v = v_cont;
                p = p_cont;
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
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                /* параметры за правой ударной волной */
                r = rr * (p_ratio + g6) / (p_ratio * g6 + 1.0);
                v = v_cont;
                p = p_cont;
            }
        }
        else {
            /* правая волна разрежения */
            shr = vr + cr;
            if (s >= shr) {
                /* параметры справа от разрыва */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                cmr = cr * pow(p_cont / pr, g1);
                str = v_cont + cmr;
                if (s <= str) {
                    /* параметры справа от контактного разрыва */
                    r = rr * pow(p_cont / pr, 1.0 / params->g);
                    v = v_cont;
                    p = p_cont;
                }
                else {
                    /* параметры внутри правой волны разрежения */
                    v = g5 * (-cr + g7 * vr + s);
                    c = g5 * (cr - g7 * (vr - s));
                    r = rr * pow(c / cr, g4);
                    p = pr * pow(c / cr, g3);
                }
            }
        }
    }
    /* формирование выходного вектора с результатом */
    r_res = r;
    v_res = v;
    p_res = p;

}

void godunov::Riman_solver(struct Parameters* params, double rl, double vl, double pl, double rr, double vr, double pr, double& p_res, double& v_res, double& r_res) {
    double cr, cl;
    double p_cont, v_cont;
    double p, v, r;
    double g = params->g;

    cl = Sound_velocity(params, pl, rl);
    cr = Sound_velocity(params, pr, rr);

    if (2.0 * (cl + cr) / (g - 1.0) <= vr - vl) {
        /* случай возникновения вакуума */
        printf("\nContact_pressure_velocity -> vacuum is generated\n");
    }
    /* итерационная процедура расчета давления и скорости газа на контактном разрыве*/
    Contact_pressure_velocity(params, pl, vl, rl, cl, pr, vr, rr, cr, p_cont, v_cont);
    /* отбор решения */
    Sample_solid_solution(params, pl, vl, rl, cl, pr, vr, rr, cr, p_cont, v_cont, 0.0, p, v, r);
    p_res = p;
    v_res = v;
    r_res = r;
}

void godunov::Diff_flux_ncons_x(struct Parameters* params, double p, double vx, double r, double& Fm, double& Fimp_x, double& Fe) {
    double m, impx, e; /* консервативные переменные */

    Convert_noncons_to_cons(params, p, vx, r, m, impx, e);
    Fm =  r * vx;
    Fimp_x = Fm * vx + p;
    Fe = (p + e) * vx;
}

void godunov::Godunov_flux_x(struct Parameters* params, double ml, double impxl, double el, double mr, double impxr, double er, double& Fm, double& Fimp_x, double& Fe) {
    double p, vx, r;
    double pl, vxl, rl;
    double pr, vxr, rr;

    Convert_cons_to_noncons(params, pl, vxl, rl, ml, impxl, el);
    Convert_cons_to_noncons(params, pr, vxr, rr, mr, impxr, er);

    /*решение задачи о распаде разрыва*/
    Riman_solver(params, rl, vxl, pl, rr, vxr, pr, p, vx, r);

    /* расчет потока Годунова по вектору неконсервативных переменных*/
    Diff_flux_ncons_x(params, p, vx, r, Fm, Fimp_x, Fe);
}

double godunov::calc_time_step(struct Parameters* params, double* x, double* m, double* impx, double* e, int time_step_number) {
    double new_step = 1000000;
    double p, vx, r, c;
    double c_step;
    double CFL = params->CFL;
    for (int i = 0; i < params->cells_number_x; ++i) {
        Convert_cons_to_noncons(params, p, vx, r, m[i], impx[i], e[i]);
        c = Sound_velocity(params, p, r);
        c_step = CFL * (x[i + 1] - x[i]) / (std::fabs(vx) + c);
        if (c_step < new_step) {
            new_step = c_step;
        }
    }
    return new_step;
}



void godunov::Init_solution_circle(struct Parameters* params, double* xc, double* p, double* vx, double* r, double* m, double* impx, double* e) {
    double p1, rho1, v1;
    double p2, rho2, v2;
    rho1 = 1.0;
    v1 = -2.0;
    p1 = 0.4;
    rho2 = 1.0;
    v2 = 2.0;
    p2 = 0.4;

    for (int i = 0; i < params->cells_number_x; ++i) {
        if (xc[i] < 0.5) {
            p[i] = p1;
            r[i] = rho1;
            vx[i] = v1;
        }
        else {
            p[i] = p2;
            r[i] = rho2;
            vx[i] = v2;
        }
        Convert_noncons_to_cons(params, p[i], vx[i], r[i], m[i], impx[i], e[i]);
    }
}


void godunov::Out(struct Parameters* params, double time, double* x, double* p, double* vx, double* r) {
    int Nx = params->cells_number_x;
    char Name_file[100];
    sprintf(Name_file, "results/%f.csv", time);
    std::ofstream SurfaceFile(Name_file, std::ios::app);
    SurfaceFile << "x;p;vx;r;e;\n";
    for (int i = 0; i < Nx; ++i) {
        SurfaceFile << x[i] << ";" << p[i] << ";" << vx[i] << ";" << r[i] << ";" << p[i] / ((params->g - 1) * r[i]) << "\n";
    }
}


void godunov::main_solver() {
    struct Parameters params;   /* структура с параметрами вычислительного эксперимента  */
    double* xc;                 /* массив x-координат центров ячеек сетки */
    double* x;                  /* массив x-координат узлов сетки */
    double* p, * vx, * r;          /* массивы неконс переменных */
    double* m, * impx, * e;        /* массивы конс переменных */
    double* m_next, * impx_next, * e_next;
    double mb, impxb, eb;        /* граничные значения конс переменных */
    double FmL, FimpxL, FeL, FmR, FimpxR, FeR; /* потоки конс переменных */

    double pl, vl, rl, pr, vr, rr;
    double ml, impxl, el, mr, impxr, er;
    double dm, dimp, de;

    int step_number = 0;
    double time = 0;
    double dt, dx;
    /* считывание файла с параметрами задачи */
    read_parameters(&params);
    int Nx = params.cells_number_x;
    int S_type = params.S_type;
    int fo = params.fo;
    int type_b_left = params.type_b_left;
    int type_b_right = params.type_b_right;
    /* выделение памяти под массивы */
    Allocate(Nx, p, vx, r, m, impx, e, m_next, impx_next, e_next);
    Allocate(Nx + 1, xc);
    Allocate(Nx + 1, x);
    /* определение координат центров ячеек сетки */
    Build_grid(&params, xc, x);
    /* начальные условия */
    Init_solution_circle(&params, xc, p, vx, r, m, impx, e);
    /* вывод в файл начального распределения*/
    Out(&params, time, xc, p, vx, r);
    /* цикл по времени*/
    while (time < params.stop_time) {
        dt = calc_time_step(&params, x, m, impx, e, step_number);
        for (int i = 0; i < Nx; ++i){
            /* ПО ОСИ X */
            /* расчет потока через левую грань ячейки */
            if (S_type == 1) {
                /* ГОДУНОВ */
                if (i != 0) {
                    ml = m[i - 1];
                    impxl = impx[i - 1];
                    el = e[i - 1];
                    mr = m[i];
                    impxr = impx[i];
                    er = e[i];
                }
                else {
                    Boundary_x(m[0], impx[0], e[0], mb, impxb, eb, params.type_b_left);
                    ml = mb;
                    impxl = impxb;
                    el = eb;
                    mr = m[i];
                    impxr = impx[i];
                    er = e[i];
                }
            }
            Godunov_flux_x(&params, ml, impxl, el, mr, impxr, er, FmL, FimpxL, FeL);

            /* расчет потока через правую грань ячейки */
            if (S_type == 1) {
                /* ГОДУНОВ */
                if (i != Nx - 1) {
                    ml = m[i];
                    impxl = impx[i];
                    el = e[i];
                    mr = m[i + 1];
                    impxr = impx[i + 1];
                    er = e[i + 1];
                }
                else {
                    Boundary_x(m[Nx - 1], impx[Nx - 1], e[Nx - 1], mb, impxb, eb, params.type_b_right);
                    ml = m[i];
                    impxl = impx[i];
                    el = e[i];
                    mr = mb;
                    impxr = impxb;
                    er = eb;
                }
            }
            Godunov_flux_x(&params, ml, impxl, el, mr, impxr, er, FmR, FimpxR, FeR);

            if (i != Nx - 1)
                dx = (xc[i + 1] - xc[i]);
            else
                dx = (xc[i] - xc[i - 1]);
            m_next[i] = m[i] - dt * (FmR - FmL) / dx;
            impx_next[i] = impx[i] - dt * (FimpxR - FimpxL) / dx;
            e_next[i] = e[i] - dt * (FeR - FeL) / dx;
            }


        for (int i = 0; i < Nx; i++) {              
            impx[i] = impx_next[i];
            m[i] = m_next[i];
            e[i] = e_next[i];
            Convert_cons_to_noncons(&params, p[i], vx[i], r[i], m[i], impx[i], e[i]);
        }
        time += dt;
        step_number++;
        /* запись в файл */
        if (step_number % fo == 0) {
       // if (time >= 0.119363) {
            Out(&params, time, xc, p, vx, r);
         //   break;
        }
    }
    /* освобождение памяти */
    Delete(x, xc, p, vx, r, m, impx, e, m_next, impx_next, e_next);
}

godunov::godunov() {
	main_solver();
}