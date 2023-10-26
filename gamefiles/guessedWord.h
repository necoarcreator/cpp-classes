#pragma once
#include <iostream>
#include <string>

using namespace std;



class GuessedWord
{

private:

    string m_truth, m_guess;

    int m_length;

    bool m_isRightGuess = false;



public:

    void setWord(string word)
    {
        this->m_truth = word;
        m_length = m_truth.size();
        m_guess += string(m_length, '-');
    }

    void print()
    {
        cout << m_truth << endl;
    }

    void getGuess(string* word)
    {
        word = &m_guess;
    }

    string getTruth()
    {
        return m_truth;
    }

    void setIsRightGuess(bool isRightGuess)
    {
        m_isRightGuess = isRightGuess;
    }

    bool getIsRightGuess()
    {
        return m_isRightGuess;
    }

    int hasLett(char buk)
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

    bool proveRight()
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

    ~GuessedWord() {}
};

