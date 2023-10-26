
#include <iostream>
#include <string>
#include <thread>
#include <chrono>
#include "dict.h"
#include "yakubovich.h"
#include "guessedWord.h"

using namespace std;



    class PlayGame
    {

    public:
        GuessedWord word;
        yakubovich LeonidArcadievich;
        dictionary communarRevolutioner;
        string new_word = communarRevolutioner.getWord();
        string answer = "n";
        string* guess;

        bool is_victory = false;

        int score;

        char buk;

        int Play()
        {
            communarRevolutioner.addWord();

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
                    }
                    score -= 1;
                    LeonidArcadievich.comment(word.getIsRightGuess());
                    word.setIsRightGuess(false);
                    buk = LeonidArcadievich.getLetter();
                }
                new_word.clear();

            
            if (not is_victory)
            {
                cout << "You lost. The right word was: ";
                cout << word.getTruth() << endl;
                cout << "\n" << endl;
                *guess = "";
                is_victory = false;
            }

            cout << "Do you want to play one more time? Enter [Y/n]: ";
                cin >> answer;
            if (answer == "Y")
            {

                Play();
            }
            else
            {
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