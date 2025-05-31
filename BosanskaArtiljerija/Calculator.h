#ifndef CALC
#define CALC

#include <list>
#include <vector>
#include <string>
#include <memory>

using namespace std;
class Calculator
{
private:
	size_t numBat;
	size_t numTarg;
	size_t numMis;

	vector<vector<string>> artyInfo;
	vector<vector<string>> targInfo;

	vector<pair<size_t, size_t>> calculatedNum; //пары (число снарядов на поражение, время на поражение в минутах)
	size_t r; //число базисных переменных 2 * numBat + numTarg
	size_t k; //число свободных переменных numTarg * numBat
	vector<vector<double>>* AMtx;
	vector<double>* BVec;
	vector<double>* CVec;
	double val;


public:
	Calculator();
	~Calculator();
	size_t countRangesAndAmount(size_t whatTarget);

	void simplexInit(double a, double b);
	void outputInfeasible(vector<vector<double> >* a, vector<double>* b, vector<double>* c);
	void reduce(vector<vector<double> >* a, vector<double>* b, vector<double>* c, double* opt, int row, int col);
	void outputUnbounded(vector<vector<double> >* a, vector<double>* b, vector<double>* c, bool minimize);
	void output(vector<vector<double> >* a, vector<double>* b, vector<double>* c, vector<double> x, double opt, int numVars);
	void reduceC(vector<vector<double> >* a, vector<double>* b, vector<double>* c, double* opt);
	void simplex(vector<vector<double>>* a, vector<double>* b, vector<double>* c, 
		vector<double> initialC, double initialOpt, int numVars, int phase, bool minimize, double _aa, double _bb);
	void gomoriInit(double a, double b);

};

#endif