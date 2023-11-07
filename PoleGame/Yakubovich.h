#pragma once

#include <iostream>
#include <vector>
#include <string>

using namespace std;

class Yakubovich
{
    int random, length_w, length_t;
    vector<string> right, wrong;
    char buk;
    string word;

public:

    Yakubovich();

    void setPhrases();

    bool checkcorr(char buk);

    char getLetter();

    void comment(bool is_right);
    ~Yakubovich();
};
