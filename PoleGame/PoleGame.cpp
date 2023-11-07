
#include <iostream>
#include <string>
#include <thread>
#include <chrono>
#include "PlayGame.h"

using namespace std;



int main()
{
    int num = 10, sum_score = 0;
    PlayGame poleChudes;

    sum_score += poleChudes.play();

    cout << "Your final score for today is " << sum_score << ". Congraduations!" << endl;


};