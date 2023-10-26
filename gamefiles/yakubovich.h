#pragma once
#include <iostream>
#include <io.h>
#include <fcntl.h>
#include <string>
#include <fstream>
#include <vector>
#include <thread>
#include <chrono>

using namespace std;


using namespace std;
class yakubovich
{
    int random, length_w = 0, length_t = 0;
    vector<string> right, wrong;
    char buk;
    string word;

public:
    void setPhrases()
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

    bool checkcorr(char buk)
    {
        if (((static_cast<int>(buk) > 64) and (static_cast<int>(buk) < 91)) or ((static_cast<int>(buk) > 96) and (static_cast<int>(buk) < 123))) //если соотв кодировке букв
        {
            return true;
        }
        else
        {
            return false;
        }

    }

    char getLetter()
    {
        cout << "Enter your guess" << endl;
        cin >> buk;
        while (not checkcorr(buk))
        {
            cout << "You've written a forbidden letter. Try again: " << endl;
            cin >> buk;
        }
        return tolower(buk);
    }

    void comment(bool is_right)
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
    ~yakubovich() {}
};