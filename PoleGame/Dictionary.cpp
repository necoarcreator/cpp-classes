#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <chrono>
#include "Dictionary.h"

using namespace std;

Dictionary::Dictionary() : random(0), length(0), word(""), answer("n") {}

void Dictionary::setDict(int num = 1)
{
    fstream buf("dict1.txt");
    if (!buf)
    {
        // то выводим следующее сообщение об ошибке и выполняем функцию exit()
        cerr << "Uh oh, dict1.txt could not be opened for reading!" << endl;
        exit(1);
    }


    while (buf)
    {
        buf >> word;
        length += 1;
        dict.push_back(word);
    }
    buf.close();
    word = "";
}

string Dictionary::getWord()
{
    setDict(1);
    srand(time(NULL));
    this_thread::sleep_for(chrono::seconds(1));
    random = rand() % length;
    return dict[random];
}

int Dictionary::addWord()
{

    cout << "Do you want to add the new word? Tap [Y/n]" << endl;
    cin >> answer;
    if (answer == "Y")
    {
        fstream bufIn("dict1.txt", ios::app);
        cout << "Insert new word: ";
        cin >> word;
        cout << word << endl;
        bufIn << "\n";
        bufIn << word;
        bufIn.close();
    }
    else
    {
        return 0;
    }
    addWord();
}

Dictionary::~Dictionary() {};
