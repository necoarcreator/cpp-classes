#include <iostream>
#include <string>
#include "GuessedWord.h"

using namespace std;

GuessedWord::GuessedWord() : m_isRightGuess(false), m_truth(""), m_guess(""), m_length(0) {}

void GuessedWord::setWord(string word)
{
    this->m_truth = word;
    m_length = m_truth.size();
    m_guess += string(m_length, '-');
}

void GuessedWord::print()
{
    cout << m_truth << endl;
}

void GuessedWord::getGuess(string* word)
{
    word = &m_guess;
}

string GuessedWord::getTruth()
{
    return m_truth;
}

void GuessedWord::setIsRightGuess(bool isRightGuess)
{
    m_isRightGuess = isRightGuess;
}

bool GuessedWord::getIsRightGuess()
{
    return m_isRightGuess;
}

int GuessedWord::hasLett(char buk)
{
    for (int i = 0; i < m_length; i++)
    {
        if (m_truth[i] == buk)
        {
            m_guess[i] = buk;
            m_isRightGuess = true;
        }


    }
    cout << m_guess << endl;
    cout << endl;

    if (m_isRightGuess)
    {
        return -1;
    }
    else
    {
        return 0;
    }

}

bool GuessedWord::proveRight()
{
    for (int i = 0; i < m_length; i++)
    {
        if (m_guess[i] == '-')
        {
            return false;
        }
    }
    if (m_truth == m_guess) {
        return true;
    }
}

void GuessedWord::clearAll()
{
    m_guess = "";
    m_truth = "";
    m_length = 0;
    return;
}


GuessedWord::~GuessedWord() {};
