#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "rectangle.h"

using namespace std;


rectangle::rectangle() : length(0), heigth(0), leftBelowPointX(0), leftBelowPointY(0)
{
}

rectangle::rectangle(float a, float b, float c, float d) : length(a), heigth(b), leftBelowPointX(c), leftBelowPointY(d)
{}


float rectangle::setProperties() {
    ifstream buf("rect.txt");
    if (!buf)
    {
        cerr << "О нет, rect.txt нельзя открыть для чтения!" << endl;
        exit(1);
    }
    while (buf)
    {
        buf >> word;
        dict.push_back(stof(word));
    }
    buf.close();
    length = dict[0];
    heigth = dict[1];
    leftBelowPointX = dict[2];
    leftBelowPointY = dict[3];
    word = "";
    cout << "gulman" << endl;
}

float* rectangle::getProperties()
{
    float* M = new float[4];
    M[0] = length;
    M[1] = heigth;
    M[2] = leftBelowPointX;
    M[3] = leftBelowPointY;
    return M;
}

rectangle::~rectangle() {
    length = 0;
    heigth = 0;
    leftBelowPointX = 0;
    leftBelowPointY = 0;
    dict.clear();
};
