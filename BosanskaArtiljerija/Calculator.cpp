#include "Calculator.h"
#include "Container.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <list>
#include <tuple>
#include <utility>
#include <exception>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <exception>
#include "ActiveArty.h"

using namespace std;

Calculator::Calculator() : numBat(0), numMis(0), numTarg(0)
{
	
	Container A;
	ifstream fileArty("arty.txt");
	try
	{
		if (fileArty.is_open())
		{
			string input;
			while (getline(fileArty, input))
			{
				stringstream sinput(input);
				char delim = ',';
				char delim2 = ';';
				string xcoord, ycoord, typearty, numinbattery, numshells, typeshells, typefuze, temperature, maxTime;
				getline(sinput, xcoord, delim);
				getline(sinput, ycoord, delim);
				getline(sinput, typearty, delim);
				getline(sinput, numinbattery, delim);
				getline(sinput, numshells, delim);
				getline(sinput, typeshells, delim);
				getline(sinput, typefuze, delim);
				getline(sinput, temperature, delim);
				getline(sinput, maxTime, delim2);
				numBat++;
				artyInfo.push_back({ xcoord, ycoord, typearty,numinbattery, numshells, typeshells, typefuze, temperature, maxTime });
				A.checkCorrectArty(typearty, static_cast<size_t>(stoi(numinbattery)), static_cast<size_t>(stoi(numshells)), 
					typeshells, typefuze, stod(temperature));
			}

		}
		else
		{
			throw runtime_error("Problem while opening file to read");
		}
		fileArty.close();
	}
	catch (exception e)
	{
		cerr << e.what() << endl;
		exit(1);
	}
	ifstream fileTarg("targets.txt");
	try
	{
		if (fileTarg.is_open())
		{
			string input;
			while (getline(fileTarg, input))
			{
				stringstream sinput(input);
				char delim = ',';
				char delim2 = ';';
				string xcoord, ycoord, typetarg, mission,entrenchedArmoured, geometry, preparing;
				getline(sinput, xcoord, delim);
				getline(sinput, ycoord, delim);
				getline(sinput, typetarg, delim);
				getline(sinput, entrenchedArmoured, delim);
				getline(sinput, mission, delim);
				getline(sinput, geometry, delim);
				getline(sinput, preparing, delim2);
				numTarg++;
				targInfo.push_back({ xcoord, ycoord, typetarg, mission, entrenchedArmoured, geometry, preparing });
				A.checkCorrectTarget(typetarg, mission, entrenchedArmoured, geometry, preparing);
			}

		}
		else
		{
			throw runtime_error("Problem while opening file to read");
		}
		fileTarg.close();
	}
	catch (exception e)
	{
		cerr << e.what() << endl;
		exit(1);
	}
	
}

Calculator::~Calculator() {
    AMtx->clear();
    delete AMtx;
    
    BVec->clear();
    delete BVec;

    CVec->clear();
    delete CVec;

}

size_t Calculator::countRangesAndAmount(size_t whatTarget)
{
    vector<ActiveArty> units;
	double xTarg = stod(targInfo[whatTarget][0]);
	double yTarg = stod(targInfo[whatTarget][1]);
	size_t numVariants = 0;
	Container A;
	size_t actualCharge = 0;
	for_each(artyInfo.begin(), artyInfo.end(), [&](auto& _n)
		{
			double xArty = stod(_n[0]);
			double yArty = stod(_n[1]);
			string name = _n[2];
			size_t numInBattery = static_cast<size_t>(stoi(_n[3]));
			size_t numShells = static_cast<size_t>(stoi(_n[4]));
			double range = pow(pow(xArty - xTarg, 2) + pow(yArty - yTarg, 2), 0.5);
			actualCharge = 1;
				//A.getCharge(name, range);
			if (range > 0 && range < A.getRange(name) && numShells > 0) {
                units.push_back({ range, stod(_n[7]), A.getType(name),A.getCaliber(name),
					numInBattery, numShells, _n[5], _n[6], name}); //создаём список активной на эту цель артиллерии

			}
		});
	//range, temperatire, <howitzer, caliber>, <numinbatt, numshellscommon>, <typeshells, typefuze>
	string targName = targInfo[whatTarget][2];
	size_t move = calculatedNum.size();
	size_t j = 0;
	vector<pair<size_t, size_t>> filling = vector<pair<size_t, size_t>>(numBat, make_pair(UINT_MAX, UINT_MAX));

	calculatedNum.insert(calculatedNum.end(), filling.begin(), filling.end()); //по дефолту активная артиллерия - вся, но для неактивной максимальные числа в парах
	switch (A.checkStationary(targName))
	{
	case (true):

		for_each(units.begin(), units.end(), [&](auto& _n)
			{
				bool isCluster = (_n.typeshells == "cluster") ? true : false;
				bool isShortened = (targInfo[whatTarget][5] == "shortened") ? true : false;
				bool isEntrenchedOrArmored = (targInfo[whatTarget][3] == "entrenched" || targInfo[whatTarget][3] == "armoured") ? true : false;
				bool isDestruction = (targInfo[whatTarget][4] == "destruction") ? true : false; 
				double temperature = _n.temperature;

				size_t actualNumShells = A.getNormStationary(_n.typesCalibers, isCluster, false,
					isShortened, isEntrenchedOrArmored, isDestruction, false, _n.range, _n.numInBattery, targName);
				size_t actualMinutes = A.measureFireRate(_n.typesCalibers, actualCharge, temperature, actualNumShells, A.getSpecificType(_n.name));
				calculatedNum[move + j] = make_pair(actualNumShells, actualMinutes);
				numVariants++;
				j++;
			});
	
		break;

	case (false):
		for_each(units.begin(), units.end(), [&](auto& _n)
			{
				bool isCluster = (_n.typeshells == "cluster") ? true : false;
				bool isEntrenchedOrArmored = (targInfo[whatTarget][3] == "entrenched" || targInfo[whatTarget][3] == "armoured") ? true : false;
				double temperature = _n.temperature;

				size_t actualNumShells = A.getNormMoving(_n.typesCalibers, isCluster, isEntrenchedOrArmored, _n.numInBattery, targName);
				size_t actualMinutes = A.measureFireRate(_n.typesCalibers, actualCharge, temperature, actualNumShells, A.getSpecificType(_n.name));
				calculatedNum[move + j] = make_pair(actualNumShells, actualMinutes);
				numVariants++;
				j++;
			});

	}

    units.clear();
	return numVariants;


}

void Calculator::simplexInit(double a, double b)
{
	vector<size_t> numVar;

	for (size_t i = 0; i < numTarg; i++)
	{
		numVar.push_back(countRangesAndAmount(i)); //рассчитываем варианты на каждую цель и добавляем их число к вектору индексов
	}
	r = 2 * numBat + numTarg; //число базисных переменных
    k = numTarg * numBat;
        //r - numTarg; //число свободных переменных


	size_t i = 0;

	AMtx = new vector<vector<double>>(r, vector<double>(k, 0));
	auto it1 = AMtx->begin();
	size_t move = 0;
	advance(it1, numBat);
	for_each(AMtx->begin(), it1, [&](auto& _n)
		{
			i = 0;
			for_each(_n.begin(), _n.end(), [&](auto& _m)
				{
					if ((move <= i) && (i < move + numTarg))
					{
						_m = calculatedNum[i].first;
					}
					else
					{
						_m = 0;
					}
					i++;
				});
			move += numTarg;
		});

	auto it2 = it1;
	advance(it2, numBat);
	move = 0;
	for_each(it1, it2, [&](auto& _n)
		{
			i = 0;
			for_each(_n.begin(), _n.end(), [&](auto& _m)
				{
					if ((move <= i) && (i < move + numTarg))
					{
						_m = calculatedNum[i].second;
					}
					else
					{
						_m = 0;
					}
					i++;
				});
			move += numTarg;
		});
	i = 0;
	for_each(it2, AMtx->end(), [&](auto& _n)
		{
			size_t j = 0;
            auto it = _n.begin();
            advance(it, numTarg * numBat);
			for_each(_n.begin(), it, [&](auto& _m)
				{
					if ((j) % numTarg - i == 0)
					{
						_m = 1;
					}
					else
					{
						_m = 0;
					}
					j++;
				});
			i++;
		});
   // for (i = 0; i < 2 * numBat; i++)
   // {
   //     (*AMtx)[i][i + numBat * numTarg] = 1.0; //добавляем базисные слагаемые
   // }

    i = 0;
	BVec = new vector<double>(r);
	auto it3 = BVec->begin();
	advance(it3, numBat);
	for_each(BVec->begin(), it3, [&](auto& _n)
		{
			_n = static_cast<double>(stof(artyInfo[i][4])); //число снарядов на батарее
			i++;
		});

	i = 0;
	auto it4 = it3;
	advance(it4, numBat);
	for_each(it3, it4, [&](auto& _n)
		{
			_n = static_cast<long long int>(stof(artyInfo[i][8])); //время до передислокации батареи
			i++;
		});

	for_each(it4, BVec->end(), [&](auto& _n)
		{
			_n = 1; //сколько батарей должно быть подключено к поражению
			i++;
		});

	i = 0;
	CVec = new vector<double>(k);
	auto it5 = CVec->begin();
	advance(it5, numBat*numTarg);
	for_each(CVec->begin(), it5, [&](auto& _n)
		{
			_n = a * calculatedNum[i].first + b * calculatedNum[i].second;
			i++;
		});

    for_each(it5, CVec->end(), [&](auto& _n)
        {
            _n = 0;
        });

	val = 0;

	simplex(AMtx, BVec, CVec, *CVec, 0, k, 1, true, a, b);

}

void Calculator::output(vector<vector<double> >* a, vector<double>* b, vector<double>* c, vector<double> x, double opt, int numVars)
{
    // open out.txt
    ofstream outfile;
    outfile.open("out.txt");
    // make sure it's open
    try
    {
        if (outfile.is_open())
        {
            // output the final simplex tableau seen
            outfile << "Final Tableau:" << endl;
            size_t i = 0;
            for_each(a->begin(), a->end(), [&](auto const& _n)
                {
                    for_each(_n.begin(), _n.end(), [&](auto const& _m)
                        {

                            // output the A matrix
                            outfile << _m << " ";
                        });
                    // output the b matrix
                    outfile << "| " << (*b)[i] << endl;
                    i++;
                });

            // output the c matrix
            for_each(c->begin(), c->end(), [&](auto const& _n)
                {
                    outfile << _n << " ";
                });
            // output the optimal value found
            outfile << "| " << opt << endl;
            outfile << endl;

            outfile << "z* = " << opt << endl;

            // output the optimal x
            outfile << "x* = (";
            for (int i = 0; i < x.size(); i++)
            {
                if (i == x.size() - 1)
                    outfile << x[i] << ")" << endl;
                else
                    outfile << x[i] << ", ";
            }

            // close out.txt
            outfile.close();
        }
        // error
        else

            throw runtime_error("Could not open output file.");
    }
    catch (exception e)
    {
        cerr << e.what() << endl;
        exit(1);
    }


}

//! Outputs an unbounded solution to "out.txt"
/*!
  \param a the A matrix
  \param b the b matrix
  \param c the c matrix
  \param minimize true if min, false if max
  \sa output(), outputInfeasible()
*/
void Calculator::outputUnbounded(vector<vector<double> >* a, vector<double>* b, vector<double>* c, bool minimize)
{
    // open out.txt
    ofstream outfile;
    outfile.open("out.txt");
    // make sure it's open
    try
    {
        if (outfile.is_open())
        {
            // output the final simplex tableau seen
            outfile << "Final Tableau:" << endl;
            size_t i = 0;
            for_each(a->begin(), a->end(), [&](auto const& _n)
                {
                    for_each(_n.begin(), _n.end(), [&](auto const& _m)
                        {

                            // output the A matrix
                            outfile << _m << " ";
                        });
                    // output the b matrix
                    outfile << "| " << (*b)[i] << endl;
                    i++;
                });
            // output the c matrix
            for_each(c->begin(), c->end(), [&](auto const& _n)
                {
                    outfile << _n << " ";
                });
            outfile << endl;

            // optimal is unbounded
            if (minimize)
                outfile << "z* = -infinity" << endl;
            else
                outfile << "z* = infinity" << endl;
            outfile << "The program is unbounded." << endl;

            // close out.txt
            outfile.close();
        }
        else throw runtime_error("Could not open output file.");
    }
    catch (exception e)
    {
        cerr << e.what() << endl;
        exit(1);
    }
        
}

//! Outputs an infeasible solution to "out.txt"
/*!
  \param a the A matrix
  \param b the b matrix
  \param c the c matrix
  \sa output(), outputUnbounded()
*/
void Calculator::outputInfeasible(vector<vector<double> >* a, vector<double>* b, vector<double>* c)
{
    // open out.txt
    ofstream outfile;
    outfile.open("out.txt");
    try
    {
        // make sure it's open
        if (outfile.is_open())
        {
            // output the final simplex tableau seen
            outfile << "Final Tableau:" << endl;
            size_t i = 0;
            for_each(a->begin(), a->end(), [&](auto const& _n)
                {
                    for_each(_n.begin(), _n.end(), [&](auto const& _m)
                        {

                            // output the A matrix
                            outfile << _m << " ";
                        });
                    // output the b matrix
                    outfile << "| " << (*b)[i] << endl;
                    i++;
                });
            // output the c matrix
            for_each(c->begin(), c->end(), [&](auto const& _n)
                {
                    outfile << _n << " ";
                });
            outfile << endl;

            // problem is infeasible
            outfile << "The program is infeasible." << endl;

            // close out.txt
            outfile.close();
        }
        // error
        else throw runtime_error("Could not open output file.");
    }
    catch (exception e)
    {
        cerr << e.what() << endl;
        exit(1);
    }
}
void Calculator::reduce(vector<vector<double> >* a, vector<double>* b, vector<double>* c, double* opt, int row, int col)
{
    // get the value needed to make the pivot element 1
    double pivdiv = 1 / (*a)[row][col];
    // reduce the pivot row with this value
    for (int i = 0; i < (*a)[0].size(); i++)
    {
        (*a)[row][i] = (*a)[row][i] * pivdiv;
    }
    (*b)[row] = (*b)[row] * pivdiv;

    // for every other row
    for (int i = 0; i < (*a).size(); i++)
    {
        if (i == row)
            continue;

        // get the value needed to make the pivot column element 0
        double coeff = (*a)[i][col] * (*a)[row][col];
        // reduce the row with this value
        for (int j = 0; j < (*a)[0].size(); j++)
        {
            (*a)[i][j] = (*a)[i][j] - ((*a)[row][j] * coeff);
        }
        (*b)[i] = (*b)[i] - ((*b)[row] * coeff);
    }
    // get the value needed to make the c element 0 in the pivot column
    double coeff = (*c)[col] * (*a)[row][col];
    // reduce the c row with this value
    for (int i = 0; i < (*c).size(); i++)
    {
        (*c)[i] = (*c)[i] - ((*a)[row][i] * coeff);
    }
    (*opt) = (*opt) - ((*b)[row] * coeff);
}

//! Performs row operations to make c elements in basic columns 0
/*!
  \param a the A matrix
  \param b the b matrix
  \param c the c matrix
  \param opt the optimal value in the current tableau
  \sa reduce()
*/
void Calculator::reduceC(vector<vector<double> >* a, vector<double>* b, vector<double>* c, double* opt)
{
    // for each row
        for (int i = 0; i < (*a).size(); i++)
        {
            // default the row operation coefficient to 1
            double coeff = 1;
            for (int j = 0; j < (*a)[0].size(); j++)
            {
                // found a 1, check for a basic column
                if ((*a)[i][j] == 1)
                {
                    bool basic = true;
                    for (int k = 0; k < (*a).size(); k++)
                    {
                        // not 0 or 1, not basic
                        if ((*a)[k][j] != 0 && k != i)
                            basic = false;
                    }
                    // basic, make sure the c element becomes 0
                    if (basic)
                        coeff = (*c)[j] * (*a)[i][j];
                }
                // reduce the c element in this column
                (*c)[j] = (*c)[j] - ((*a)[i][j] * coeff);
            }
            // apply the reduction to the optimal value
            (*opt) = (*opt) - ((*b)[i] * coeff);
        }
        return;
}

//! Performs the Simplex algorithm on a given tableau
/*!
  \param a the A matrix
  \param b the b matrix
  \param c the c matrix
  \param initialC the original c matrix (for calculating optimal)
  \param initialOpt the optimal value to begin the algorithm with
  \param numVars the number of variables in the program
  \param phase 1 or 2, corresponding to what phase of the Simplex Method we are on
  \param minimize true if min, false if max
*/
void Calculator::simplex(vector<vector<double>>* a, vector<double>* b, vector<double>* c,
    vector<double> initialC, double initialOpt, int numVars, int phase, bool minimize, double _aa, double _bb)
{
    // Phase 1
    if (phase == 1)
    {
        // make the phase 1 c row
        vector<double> p1c;
        // make original variables 0
        for (int i = 0; i < (*a)[0].size(); i++)
        {
            p1c.push_back(0.0);
        }
        

        // add as many 1s as rows

            for (int i = 0; i < (*a).size(); i++)
            {
                if (i < (*a).size() - numTarg) p1c.push_back(1.0);
                else p1c.push_back(-1.0);
                // add an identity matrix to A

                for (int j = 0; j < (*a).size(); j++)
                {
                    if (j < (*a).size() - numTarg)
                    {
                        if (j == i)
                            (*a)[i].push_back(1.0);
                        else
                            (*a)[i].push_back(0.0);
                    }
                    else
                    {
                        if (j == i)
                            (*a)[i].push_back(-1.0);
                        else
                            (*a)[i].push_back(0.0);
                    }
                }
                
            }
           
        double opt = 0.0;
        // reduce the new c row so basic variables have 0 cost
        reduceC(a, b, &p1c, &opt);

        // run the Simplex algorithm
        bool stop = false;
        while (!stop)
        {
            // get the pivot column
            double min_r = 0;
            int piv_col = 0;
            for (int i = 0; i < p1c.size(); i++)
            {
                // check for the smallest c < 0
                if (p1c[i] < min_r)
                {
                    min_r = p1c[i];
                    piv_col = i;
                }
            }

            // no c < 0, we're done
            if (min_r >= 0)
            {
                stop = true;
                continue;
            }

            // get the pivot row
            vector<double> ratios;
            // for each row
            for (int i = 0; i < (*a).size(); i++)
            {
                if ((*a)[i][piv_col] == 0)
                    ratios.push_back(-1);
                else
                {
                    // denegerate loop ahead, don't pivot here
                    if ((*a)[i][piv_col] <= 0 && (*b)[i] == 0)
                        ratios.push_back(-1);
                    else
                        ratios.push_back((*b)[i] / (*a)[i][piv_col]);
                }
            }
            double min_ratio = ratios[0];
            int piv_row = 0;
            for (int i = 1; i < ratios.size(); i++)
            {
                if (ratios[i] >= 0 && (ratios[i] < min_ratio || min_ratio < 0))
                {
                    min_ratio = ratios[i];
                    piv_row = i;
                }
            }

            // no positive ratios, stop here
            if (min_ratio < 0)
            {
                stop = true;
                continue;
            }

            // reduce the tableau on the pivot row and column
            reduce(a, b, &p1c, &opt, piv_row, piv_col);

            // return to step 1 of the Simplex algorithm
            long long int indFrac = -1;
            double maxFracPart = 0.0;
            bool isInteger = true;
            for (auto x_i : (*b))
            {
                if (x_i != floor(x_i))
                {
                    isInteger = false;
                    break;
                }
            }

            if (!isInteger)
            {
                for (size_t i = 0; i < (*b).size(); i++)
                {
                    double fracPart = 0.;
                    if ((*b)[i] > 0) fracPart = (*b)[i] - floor((*b)[i]);
                    else fracPart = (*b)[i] - ceil((*b)[i]);

                    if (fracPart > maxFracPart) {
                        maxFracPart = fracPart;
                        indFrac = i;
                    }
                }

                r++;
                k++;
                numVars++;
                AMtx->push_back(vector<double>((*AMtx)[0].size() + 1, 0.));
                for (size_t i = 0; i < r - 1; i++)
                {
                    (*AMtx)[i].push_back(0.);
                }
                //ещё ограничение и переменная


                for (size_t j = 0; j < (*AMtx)[0].size() - 1; j++)
                {
                    if ((*AMtx)[indFrac][j] > 0) (*AMtx)[r - 1][j] = (*AMtx)[indFrac][j] - floor((*AMtx)[indFrac][j]);
                    else (*AMtx)[r - 1][j] = (*AMtx)[indFrac][j] - ceil((*AMtx)[indFrac][j]);
                }
                double bVal = 0;
                if ((*b)[indFrac] > 0) bVal = (*b)[indFrac] - floor((*b)[indFrac]);
                else  bVal = (*b)[indFrac] - ceil((*b)[indFrac]);

                BVec->push_back(abs(bVal));
                if (bVal > 0)
                {
                    (*AMtx)[r - 1][(*AMtx)[0].size() - 1] = 1.;
                }
                else
                {
                    (*AMtx)[r - 1][(*AMtx)[0].size() - 1] = -1.;
                }

                CVec->push_back(0.);

                
                for (size_t j = 0; j < k; j++)
                {
                    if (((*CVec)[j] > 0.) && ((*CVec)[j] < 1e-5)) (*CVec)[j] = 0.;
                }

                for (size_t i = 0; i < r; i++)
                {
                    if (((*BVec)[i] > 0) && ((*BVec)[i] < 1e-5)) (*BVec)[i] = 0.;
                    for (size_t j = 0; j < (*AMtx)[0].size(); j++)
                    {
                        if (((*AMtx)[i][j] > 0) && ((*AMtx)[i][j] < 1e-5)) (*AMtx)[i][j] = 0.;

                    }
                }
            }
            
        }

        // phase 1 done, cut off the added variables
        for (int i = 0; i < (*a).size(); i++)
        {
            (*a)[i].resize(numVars);
        }
        for (size_t j = initialC.size(); j < CVec->size(); j++)
        {
            initialC.push_back((*CVec)[j]);
        }
        // reset optimal value
        opt = 0;
        // reduce the new c row so basic variables have 0 cost
        reduceC(a, b, c, &opt);
        // go to phase 2
        simplex(a, b, c, initialC, opt, numVars, 2, minimize, _aa, _bb);
    }
    // Phase 2
    else
    {
        double opt = initialOpt;
        
        // run the Simplex algorithm
        bool stop = false;
        while (!stop)
        {
            // get the pivot column
            double min_r = 0;
            int piv_col = 0;
            for (int i = 0; i < (*c).size(); i++)
            {
                // check for the smallest c < 0
                if ((*c)[i] < min_r)
                {
                    min_r = (*c)[i];
                    piv_col = i;
                }
            }

            // no c < 0, we're done
            if (min_r >= 0)
            {               
                stop = true;
                // check for infeasibility
                bool infeasible = true;
                for (int i = 0; i < (*c).size(); i++)
                {
                    // any c > 0, feasible
                    if ((*c)[i] > 0)
                    {
                        infeasible = false;
                        break;
                    }
                }

                if (infeasible)
                {
                    outputInfeasible(a, b, c);
                    continue;
                }

                // get x* from the matrix
                vector<double> x;
                // for each column
                for (int i = 0; i < numVars; i++)
                {
                    // check if x_i is basic
                    bool basic = true;
                    int basic_i = 0;
                    for (int j = 0; j < (*a).size(); j++)
                    {
                        // not 0 or 1, not basic
                        if ((*a)[j][i] != 0 && (*a)[j][i] != 1)
                        {
                            basic = false;
                        }
                        // get the basic index
                        if ((*a)[j][i] == 1)
                            basic_i = j;
                    }

                    // basic x, get the b value
                    if (basic)
                        x.push_back((*b)[basic_i]);
                    // nonbasic, 0
                    else
                        x.push_back(0.0);
                }
            
               
                // get optimal
                opt = 0;
                // multiply x values with cost function, sum up
                for (int i = 0; i < x.size(); i++) {
                    opt += x[i] * initialC[i];
                }

                

                // min function, just output
                if (minimize)
                {
                    output(a, b, c, x, opt, numVars);
                    
                }

                // max function, output with negative optimal
                else if (!minimize)
                {
                    output(a, b, c, x, -1 * opt, numVars);
                    
                }
                continue;
            }

            // get the pivot row
            vector<double> ratios;
            // for each row
            for (int i = 0; i < (*a).size(); i++)
            {
                if ((*a)[i][piv_col] == 0)
                    ratios.push_back(-1);
                else
                {
                    // denegerate loop ahead, don't pivot here
                    if ((*a)[i][piv_col] <= 0 && (*b)[i] == 0)
                        ratios.push_back(-1);
                    else
                        ratios.push_back((*b)[i] / (*a)[i][piv_col]);
                }
            }
            double min_ratio = ratios[0];
            int piv_row = 0;
            for (int i = 1; i < ratios.size(); i++)
            {
                if (ratios[i] >= 0 && (ratios[i] < min_ratio || min_ratio < 0))
                {
                    min_ratio = ratios[i];
                    piv_row = i;
                }
            }

           

            // no positive ratios, problem is unbounded
            if (min_ratio < 0)
            {
                stop = true;
                outputUnbounded(a, b, c, minimize);
                continue;
            }

            // reduce the tableau on the pivot row and column
            reduce(a, b, c, &opt, piv_row, piv_col);

            // return to step 1 of the Simplex algorithm
        }
    }
}

void Calculator::gomoriInit(double a, double b)
{

}