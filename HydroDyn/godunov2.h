#pragma once

#include "krootoy.h"

#include "tests.h"

class godunov: public krootoy {
	public:
	
	Grid* imp;  // импульс
	Grid* imp_new; // импульс на след шаге
	
	Grid* p_an;  // давление для аналитики задачи Сода
	Grid* v_an; // скорость
	Grid* rho_an; // плотность
	Grid* m_an; // масса
	Grid* imp_an; // импульс
	Grid* e_an; // поолная энергия


	godunov(string file_path);
	~godunov();
	void Convert_cons_to_noncons(double& p, double& v, double& rho, double m, double imp, double e);
	void Convert_cons_to_noncons(Grid* p, Grid* v, Grid* rho, Grid* m, Grid* imp, Grid* e);
	void Convert_noncons_to_cons(double p, double v, double rho, double& m, double& imp, double& e);
	void Convert_noncons_to_cons(Grid* p, Grid* v, Grid* rho, Grid* m, Grid* imp, Grid* e);

	double Get_dt ();
	double Sound_velocity (double p_sound_velocity, double rho_sound_velocity);
	
	void Solver_godunov(int test = 1, string file_path_res = "\\results", string file_path_an = "\\results_new", string format = ".csv", char _delim = '/t', int fo = 100);
	
	void Godunov_flux(double ml, double impxl, double el, double mr, double impxr, double er, double& Fm, double& Fimp_x, double& Fe);
	void Diff_flux_ncons_x(double p_flux, double vx_flux, double r, double& Fm, double& Fimp_x, double& Fe);
	void Riman_solver(double rl, double vl, double pl, double rr, double vr, double pr, double& p_res, double& v_res, double& r_res);
	void Sample_solid_solution(double pl, double vl, double rl, double cl, double pr, double vr, double rr, double cr,
		double p_cont, double v_cont, double s, double& p_res, double& v_res, double& r_res);
	void Contact_pressure_velocity(double pl, double vl, double rl, double cl,
		double pr, double vr, double rr, double cr, double& p_cont, double& v_cont);
	double Pressure_initial_guess(double pl, double vl, double rl,
		double cl, double pr, double vr, double rr, double cr);
	void Calc_F_and_DF(double curr_press, double p_F_and_DF, double v_F_and_DF, double r, double c_F_and_DF, double* F, double* DF);
	void BoundaryWallNonCons(bool type, Grid& p_boundary, Grid& v_boundary, Grid& rho_boundary);
	void BoundaryWallCons(bool type, Grid& m_boundary, Grid& imp_boundary, Grid& e_boundary);
	void BoundaryOpenNonCons(bool type, Grid& p_boundary, Grid& v_boundary, Grid& rho_boundary);
	void BoundaryOpenCons(bool type, Grid& m_boundary, Grid& imp_boundary, Grid& e_boundary);
	void writeResults(string file_path_res, string file_path_an, string format, char _delim, double _time, int test); //переопределяем, т.к. требуется иная запись
	void Analytic(double time, double l, double x, double ml, double impl, double el, double mr, double impr, double er,
					double& p_out, double& v_out, double& rho_out, double& m_out, double& imp_out, double& e_out);


	void calc_ei();

	void Kolgan_right(double qi, double qi_1, double qi1, double qi2, double xi, double xi_1, double xi1, double xi2, double& ql, double& qr, int mod_type);
	void Kolgan_left(double qi, double qi_1, double qi1, double qi2, double xi, double xi_1, double xi1, double xi2, double& ql, double& qr, int mod_type);
	void Solver_godunov_kolgan(int test = 1, string file_path_res = "\\results", string file_path_an = "\\results_new", string format= ".csv", char _delim = '/t', int fo = 100, int mod_type = 1);
	double min_mod_1984 (double a, double b);
	double min_mod_1975 (double a, double b);
	double min_mod_1972 (double a, double b);
	void Solver_godunov_WENO(int test = 1, string file_path_res = "\\results", string file_path_an = "\\results_new", string format = ".csv", char _delim = '/t', int fo = 100);
	
	void WENO_reconstructRight(Grid v, Grid& Flow);
	void WENO_reconstructLeft(Grid v, Grid& Flow);
	void WENO_flux(double p, double v, double rho, double c_max, double Fm, double Fimp, double Fe, 
																								double& Fm_minus, double& Fimp_minus, double& Fe_minus, double& Fm_plus, double& Fimp_plus, double& Fe_plus);

	 

	void Solver_mccormak(int test, string file_path_res, string file_path_an, string format, char _delim, int fo);
	void Mccormak_flux(vector<Grid> u, vector<Grid>& Flux, size_t i, bool leftOrRight);
	void Diffusion(vector<Grid> u, vector<Grid>& Du);
	void Antidiffusion(vector<Grid>& Nu);
};