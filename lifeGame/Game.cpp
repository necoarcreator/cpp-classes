#include "Game.h"
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <random>

#include "Game.h"
#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <random>
#include <algorithm>

Game::Game(std::string _input, size_t _sizeField,
    size_t _numToDieAlone, size_t _numToDieTogether, size_t _numBirth)
    : sizeField(_sizeField),
    numToDieAlone(_numToDieAlone),
    numToDieTogether(_numToDieTogether),
    numBirth(_numBirth),
    wereMoves(false) {

    field = std::vector<std::vector<Cell>>(sizeField, std::vector<Cell>(sizeField));

    for (size_t i = 0; i < sizeField; ++i) {
        for (size_t j = 0; j < sizeField; ++j) {
            if (_input[i * sizeField + j] == '1') {
                field[i][j] = Cell(true);
                aliveCells.push_back(std::make_pair(i, j));
            }
            else {
                field[i][j] = Cell(false);
            }
        }
    }
}

Game::~Game() {
    std::cout << "Game ended!" << std::endl;
}

size_t Game::returnAlive(std::pair<size_t, size_t> coord) {
    size_t count = 0;
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            int x = coord.first + dx;
            int y = coord.second + dy;
            if (x >= 0 && x < (int)sizeField && y >= 0 && y < (int)sizeField) {
                if (!(dx == 0 && dy == 0)) { // не считаем саму себя
                    if (field[x][y].state)
                        ++count;
                }
            }
        }
    }
    return count;
}

void Game::makeStep() {
    std::vector<std::pair<size_t, size_t>> newBorn;
    std::vector<std::pair<size_t, size_t>> toDie;

    for (auto it = aliveCells.begin(); it != aliveCells.end(); ++it) {
        auto coord = *it;
        size_t neighbors = returnAlive(coord);

        if (neighbors < numToDieAlone || neighbors > numToDieTogether) {
            toDie.push_back(coord);
        }
        else {
            field[coord.first][coord.second].nextState = true;

            if (neighbors == numBirth) {
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_int_distribution<> dis(-1, 1);
                int dx = dis(gen);
                int dy = dis(gen);

                while (dx == 0 && dy == 0) {
                    dx = dis(gen);
                    dy = dis(gen);
                }

                int nx = coord.first + dx;
                int ny = coord.second + dy;

                if (nx >= 0 && nx < (int)sizeField && ny >= 0 && ny < (int)sizeField) {
                    if (!field[nx][ny].state) {
                        field[nx][ny].nextState = true;
                        newBorn.push_back(std::make_pair(nx, ny));
                    }
                }
            }
        }
    }

    // Применяем смерть
    for (auto& coord : toDie) {
        field[coord.first][coord.second].nextState = false;
        aliveCells.erase(std::remove(aliveCells.begin(), aliveCells.end(), coord), aliveCells.end());
    }

    // Добавляем новые клетки
    for (auto& coord : newBorn) {
        field[coord.first][coord.second].nextState = true;
        aliveCells.push_back(coord);
    }
}

void Game::refreshCells() {
    wereMoves = false;

    for (size_t i = 0; i < sizeField; ++i) {
        for (size_t j = 0; j < sizeField; ++j) {
            if (field[i][j].state != field[i][j].nextState) {
                wereMoves = true;
            }
            field[i][j].state = field[i][j].nextState;
        }
    }
}

void Game::printField() {
    for (const auto& row : field) {
        for (const auto& cell : row) {
            std::cout << (cell.state ? " @ " : " . ");
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

size_t Game::getAliveNum() {
    return aliveCells.size();
}

bool Game::getWereMoves() {
    return wereMoves;
}


