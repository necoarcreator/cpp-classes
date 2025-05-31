#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <fstream>
//#include "krootoy.cpp"
/* ПАРАМЕТРЫ */
struct PVR{
    double p;
    double v;
    double r;
};

/* ВХОДНЫЕ ДАННЫЕ */
struct godunov {
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
	
	
	void read_godunov(struct Parameters* params);
	void Convert_cons_to_noncons(struct Parameters* params, double& p, double& vx, double& r, double& m, double& impx, double& e);
	void Convert_noncons_to_cons(struct Parameters* params, double& p, double& vx, double& r, double& m, double& impx, double& e);
	double Sound_velocity(struct Parameters* params, double p, double r);
	void Build_grid(struct Parameters* params, double* xc, double* x);
	void Boundary_x(double p, double vx, double r, double& pb, double& vxb, double& rb, int b_type);
	void Calc_F_and_DF(struct Parameters* params, double curr_press, double p, double v, double r, double c, double* F, double* DF);
	double Pressure_initial_guess(struct Parameters* params, double pl, double vl, double rl,
	double cl, double pr, double vr, double rr, double cr);
	void Contact_pressure_velocity(struct Parameters* params, double pl, double vl, double rl, double cl,
    double pr, double vr, double rr, double cr, double& p_cont, double& v_cont);
	void Sample_solid_solution(struct Parameters* params, double pl, double vl, double rl, double cl, double pr, double vr, double rr, double cr,
    double p_cont, double v_cont, double s, double& p_res, double& v_res, double& r_res);
	void Riman_solver(struct Parameters* params, double rl, double vl, double pl, double rr, double vr, double pr, double& p_res, double& v_res, double& r_res);
	void Diff_flux_ncons_x(struct Parameters* params, double p, double vx, double r, double& Fm, double& Fimp_x, double& Fe);
	void Godunov_flux_x(struct Parameters* params, double ml, double impxl, double el, double mr, double impxr, double er, double& Fm, double& Fimp_x, double& Fe);
	double calc_time_step(struct Parameters* params, double* x, double* m, double* impx, double* e, int time_step_number);
	void Init_solution_circle(struct Parameters* params, double* xc, double* p, double* vx, double* r, double* m, double* impx, double* e);
	void Out(struct Parameters* params, double time, double* x, double* p, double* vx, double* r);
	void main_solver();
	godunov();
};

