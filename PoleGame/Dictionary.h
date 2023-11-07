#pragma once

#include <iostream>
#include <vector>
#include <string>

using namespace std;

class Dictionary
{
    int random, length;
    vector<string> dict;
    string word;

public:
    Dictionary();
    string answer;

    void setDict(int num);

    string getWord();
    int addWord();
    ~Dictionary();
};
