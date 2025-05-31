#include "Game.h"
#include <string>
#include <thread>
using namespace std;

//Слепов А., Б21-221

int main()
{
    string s = "0000000000000000000000000000000000001110000000000000000000000000";

    Game A(s, 10, 2, 4, 3);

    size_t numAliveNow = 1;
    size_t numAlivePrev = 0;
    bool wereMovesHere = true;
    while (numAliveNow != 0)
    {
        numAlivePrev = numAliveNow;
        A.makeStep();
       
        A.printField();
        std::this_thread::sleep_for(std::chrono::seconds(1));
        numAliveNow = A.getAliveNum();

        A.refreshCells();

        wereMovesHere = A.getWereMoves();
    }




    return 0;
}
