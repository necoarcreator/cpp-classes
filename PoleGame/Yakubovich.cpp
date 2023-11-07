#include <iostream>
#include <io.h>
#include <fcntl.h>
#include <string>
#include <fstream>
#include <vector>
#include <thread>
#include <chrono>
#include "Yakubovich.h"

using namespace std;

Yakubovich::Yakubovich() : random(0), length_w(0), length_t(0), buk('0'), word("") 
    {
    setPhrases();
    }


void Yakubovich::setPhrases()
{
    fstream tr("commright.txt");
    fstream fls("commwrong.txt");
    if (!tr)
    {
        // то выводим следующее сообщение об ошибке и выполняем функцию exit()
        cerr << "Uh oh, commright.txt could not be opened for reading!" << endl;
        exit(1);
    }
    if (!fls)
    {
        // то выводим следующее сообщение об ошибке и выполняем функцию exit()
        cerr << "Uh oh, commwrong.txt could not be opened for reading!" << endl;
        exit(1);
    }
    while (tr)
    {
        getline(tr, word);
        length_t += 1;
        right.push_back(word);
    }
    tr.close();
    word = "";

    while (fls)
    {
        getline(fls, word);
        length_w += 1;
        wrong.push_back(word);
    }
    fls.close();
    word = "";

}
char Yakubovich::getLetter()
{
    cout << "Enter your guess" << endl;
    string inpt;
    cin >> inpt;
    if (inpt.size() > 1)
    {
        cout << "Insert only one letter!" << endl;
        getLetter();
    }

    buk = inpt[0];
    inpt.clear();

    int code = static_cast<int>(buk);

    if ((code < 65) || (code > 90) && (code < 97) || (code > 122))
    {
        cout << "Don't insert non-alphabetic symbols!" << endl;
        getLetter();
    }

    if (code < 97)
    {
        return tolower(buk);
    }

    return buk;
}

void Yakubovich::comment(bool is_right)
{
    setPhrases();
    if (is_right)
    {
        srand(time(NULL));
        this_thread::sleep_for(chrono::seconds(1));
        random = rand() % length_t;
        cout << right[random] << endl;
    }
    else
    {
        srand(time(NULL));
        this_thread::sleep_for(chrono::seconds(1));
        random = rand() % length_w;
        cout << wrong[random] << endl;

    }
}
Yakubovich::~Yakubovich() {};

