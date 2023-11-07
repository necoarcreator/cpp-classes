#pragma once

#include <iostream>
#include <string>
#include <thread>
#include <chrono>

#include "Dictionary.h"
#include "Yakubovich.h"
#include "GuessedWord.h"

using namespace std;



class PlayGame
{
    GuessedWord *wordd;
    Yakubovich *LeonidArcadievich;
    Dictionary *communarRevolutioner;
    string new_word;
    string answer;
    string* guess;
    int count;

    public:

    PlayGame();

    bool is_victory;

    int score;

    char bukv;

    int play();
};

//поменять на большие букв, & вместо самих объектов