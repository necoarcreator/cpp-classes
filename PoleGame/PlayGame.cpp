
#include <iostream>
#include <string>
#include <thread>
#include <chrono>
#include "PlayGame.h"

using namespace std;


PlayGame::PlayGame() :  is_victory(false), bukv('0'), score(0), count(0)

{
    communarRevolutioner = new Dictionary();
    LeonidArcadievich = new Yakubovich();
    wordd = new GuessedWord();

};

int PlayGame::play()
{
    communarRevolutioner->addWord();

    new_word = communarRevolutioner->getWord();

    wordd->setWord(new_word);

    bukv = LeonidArcadievich->getLetter();

    score = new_word.size();

    for (int i = 0; i < new_word.size(); i++)
    {
        count = i;
        i += wordd->hasLett(bukv);
        if (i - count)
        {
            score += i - count;
        }

        wordd->getGuess(guess);
        is_victory = wordd->proveRight();
        if (is_victory)
        {
            cout << "You won. The right word was: ";
            cout << wordd->getTruth() << endl;
            cout << "\n" << endl;
            break;
        }
        LeonidArcadievich->comment(wordd->getIsRightGuess());
        wordd->setIsRightGuess(false);
        bukv = LeonidArcadievich->getLetter();
    }
    new_word.clear();


    if (not is_victory)
    {
        cout << "You lost. The right word was: ";
        cout << wordd->getTruth() << endl;
        cout << "\n" << endl;
        is_victory = false;
    }

    cout << "Do you want to play one more time? Enter [Y/n]: ";
    cin >> answer;
    if (answer == "Y")
    {
        is_victory = false;
        wordd->clearAll();

        play();
    }
    else
    {
        return score;
    }
}
