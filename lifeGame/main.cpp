#include "Game.h"
#include <thread>
#include <chrono>

int main() {
    // Вводное состояние: крестовина из живых клеток
    std::string initialState =
        "00000"
        "00100"
        "01110"
        "00100"
        "00000";

    Game game(initialState, 5, 2, 3, 3); // numToDieAlone=2, numToDieTogether=3, numBirth=3

    for (int step = 0; step < 20; ++step) {
        system("clear"); // или system("cls") на Windows
        game.printField();
        game.makeStep();
        game.refreshCells();

        if (!game.getWereMoves()) {
            std::cout << "No more moves!" << std::endl;
            break;
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(300));
    }

    return 0;
}