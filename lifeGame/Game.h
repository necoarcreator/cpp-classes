#ifndef _GAME_
#define _GAME_

#include "cell.h"
#include <iostream>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <functional>
using namespace std;

class Game {
private:
    std::vector<std::vector<Cell>> field;
    std::vector<std::pair<size_t, size_t>> aliveCells;
    size_t sizeField;
    size_t numToDieAlone;
    size_t numToDieTogether;
    size_t numBirth;
    bool wereMoves;

public:
    Game(std::string input, size_t _sizeField,
        size_t _numToDieAlone, size_t _numToDieTogether, size_t _numBirth);

    ~Game();
    size_t returnAlive(std::pair<size_t, size_t> coord);
    void makeStep();
    void printField();
    size_t getAliveNum();
    void refreshCells();
    bool getWereMoves();
};

#endif