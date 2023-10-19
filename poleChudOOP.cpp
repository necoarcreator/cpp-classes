
#include <iostream>
#include <string>
#include <thread>
#include <chrono>

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


};
class Yakubovich
{
public:

    int random;

    char buk;

     string rightLett[10] = { "Absolutely!", "As doctor recommended!", "Wonderful!",
     "Straight facts.", "Delightful!", "You did it!", "Right guess.", "I knew you're a master in it.",
     "Right!", "Truth." };

         string wrongLett[10] = { "Hilarious!", "Better luck next time.", "Never been in there!",
             "You're lying to me!", "Wrong!", "Forget about it.", "No,no,no,no,no!", "Do you really think it is so?",
             "No way.", "Absolutely not!" };



    char getLetter()
    {
        cout << "Enter your guess" << endl;
        cin >> buk;
        return buk;
    }

      void comment(bool is_right)
      {
          if (is_right)
          {
              srand(time(NULL));
              this_thread::sleep_for(chrono::seconds(1));
              random = rand() % 10;
              cout << rightLett[random] << endl;
         }
          else
          {
             srand(time(NULL));
              this_thread::sleep_for(chrono::seconds(1));
              random = rand() % 10;
              cout << wrongLett[random] << endl;

          }
      }

};
class Dictionary
{
    int random;
    string dict[20] = { "Kurdistan", "Anarchy", "communism", "exploitation", "expropriation",
        "Marx", "Engels", "Lenin", "Stalin", "Mao", "Comandante", "Subcomandante",
        "Revolution", "Commune", "Aurora", "Mayakovski", "Trotskii", "Independence",
        "Glory", "Freedom" };
    int lth = 20;

public:

    string getWord()
    {
        srand(time(NULL));
        this_thread::sleep_for(chrono::seconds(1));
        random = rand() % lth;
        return dict[random];
    }

};

    class PlayGame
    {

    public:
        GuessedWord word;
        Yakubovich LeonidArcadievich;
        Dictionary communarRevolutioner;
        string new_word = communarRevolutioner.getWord();

        string* guess;

        bool is_victory = false;

        int score;

        char buk;

        int Play()
        {
            word.setWord(new_word);

            buk = LeonidArcadievich.getLetter();

            score = new_word.size();

                for (int i = 0; i < new_word.size(); i++)
                {
                    i += word.hasLett(buk);
                    word.getGuess(guess);
                    is_victory = word.proveRight();
                    if (is_victory)
                    {
                        cout << "You won. The right word was: ";
                        cout << word.getTruth() << endl;
                        cout << "\n" << endl;
                        *guess = "";
                        is_victory = false;
                        break;
                        return score;
                    }
                    score -= 1;
                    LeonidArcadievich.comment(word.getIsRightGuess());
                    word.setIsRightGuess(false);
                    buk = LeonidArcadievich.getLetter();
                    
                    
                }

            
            if (not is_victory)
            {
                cout << "You lost. The right word was: ";
                cout << word.getTruth() << endl;
                cout << "\n" << endl;
                *guess = "";
                is_victory = false;
                return score;
            }
        }
    };




    int main()
    {
        int num = 10, sum_score;
        PlayGame poleChudes;

        for (int i = 0; i < num; i++)
        {
           sum_score += poleChudes.Play();
        }

        cout << "Your final score for today is" << sum_score << ". Congraduations!" << endl;
    };