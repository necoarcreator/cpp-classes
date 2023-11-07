#pragma once

#include <iostream>
#include <string>

using namespace std;



class GuessedWord
{

private:

    string m_truth, m_guess;

    int m_length;

    bool m_isRightGuess;



public:
    GuessedWord();
    void setWord(string word);

    void print();

    void getGuess(string* word);

    string getTruth();

    void setIsRightGuess(bool isRightGuess);

    bool getIsRightGuess();

    int hasLett(char buk);

    bool proveRight();
    
    void clearAll();

    ~GuessedWord();
};
