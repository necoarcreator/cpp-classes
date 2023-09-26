
#include <iostream>
#include <string>

using std::string;
using std::cin;
using std::cout;
using std::endl;

class GuessedWord
{

string m_truth, m_guess;

int m_tries, i, j;

char m_buk;

bool m_isRightGuess = false;


public:

    void setWord(string word)
    {
        m_truth = word;
        m_tries = m_truth.size();
        m_guess += string(m_tries, '-');
    }

    void lose()
    {   
        cout << "You lost. The right word was: ";
        cout << m_truth << endl;
        cout << "\n" << endl;
    }
    
    void win()
    {
        cout << "You won. The right word was: ";
        cout << m_truth << endl;
        cout << "\n" << endl;
    }

    void printGuess()
    {
        cout << m_guess << endl;
        cout << "\n" << endl;
    }

    void playGame(char buk)
    {
        cout << "Guess a letter!" << endl;
        printGuess();
        hasLett(buk);
        m_isRightGuess = false;
        printGuess();
        char bukNew;

        if (m_tries == 0)
        {
            lose();
        }

        else if (m_truth == m_guess)
        {
            win();
        }
        else
        {
            cin >> bukNew;
            playGame(bukNew);
        }

     }



    void hasLett(int buk)
    {
        m_buk = buk;
        for (i = 0; i < m_truth.size(); i++)
        {
            if (m_truth[i] == m_buk)
            {
                m_guess[i] = m_buk;
                m_isRightGuess = true;
            }
            if (not m_isRightGuess)
            {
                m_tries -= 1;
            }
        }
    }

};



int main()
{
    string word;

    char buk;
    
    GuessedWord slovo;

    cout << "Write a word" << endl;

    cin >> word;

    slovo.setWord(word);

    cout << "The game has started. Guess a letter!" << endl;
    cout << "\n" << endl;

    slovo.printGuess();

    cin >> buk;

    slovo.playGame(buk);

}
