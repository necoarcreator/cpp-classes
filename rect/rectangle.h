#pragma once

using namespace std;


class rectangle
{private:
    float length, heigth;
    float leftBelowPointX, leftBelowPointY;
    string word;
    vector<float> dict;
public:
    rectangle();
    rectangle(float a, float b, float c, float d);

    float setProperties();
    float* getProperties();

    ~rectangle();

};
