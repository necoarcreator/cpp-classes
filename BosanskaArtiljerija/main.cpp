
#include <iostream>
#include "Calculator.h"
#include <QApplication>
using namespace std;


int main(int argc, char** argv)
{
	double a = 1.0, b = 10.0;
	Calculator A;
	A.simplexInit(a, b);


	QApplication app(argc, argv);
	return app.exec();

	return 0;
}
