#include <iostream>
#include <conio.h>
#include <windows.h>
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
        // то выводим следующее сообщение об ошибке и выполн€ем функцию exit()
        cerr << "Uh oh, commright.txt could not be opened for reading!" << endl;
        exit(1);
    }
    if (!fls)
    {
        // то выводим следующее сообщение об ошибке и выполн€ем функцию exit()
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
    SetConsoleOutputCP(1251);
    cout << "Enter your guess" << endl;
    wstring inpt;
    wcin >> inpt;
    if (inpt.size() > 1)
    {
        cout << "Insert only one letter!" << endl;
        getLetter();
    }
    buk = inpt[0];

    int code = static_cast<int>(buk);

    try {
        string kyr = "јЅ¬√ƒ≈®∆«»… ЋћЌќѕ–—“”‘’÷„ЎўЏџ№ЁёяабвгдеЄжзийклмнопрстуфхцчшщъыьэю€";
        for (auto a : kyr)
        {
            if (code == static_cast<int>(a))
            {
                throw 1; //kyrillic letters
            }
        }
        if (inpt.length() > 2)
        {
            throw 2; //long string
        }
        else if (code == 32) //space symbol
        {
            throw 3;
        }
        else if ((code > 47) && (code < 58))
        {
            throw 4; //numbers
        }
        else if (((code < 65) || (code > 90)) && ((code < 97) || (code > 122)))
        {
            throw "Don't insert non-alphabetic symbols!";    //non-alpha symbols        
        }
        else
        {
            if (((code > 64) || (code < 91)) && ((code > 96) || (code < 123)))
            {
            }
            else
            {
                throw 1.0;
            }
        }
    }
    catch (const int exeption_num)
    {
        if (exeption_num == 1)
        {
            cout << "Error code:" << exeption_num << ". You shouldn't kirillic letters" << endl;
        }
        else if (exeption_num == 2)
        {
            cout << "Error code:" << exeption_num << ". You shouldn't write more than one symbol" << endl;
        }
    
        else if (exeption_num == 3)
        {
            cout << "Error code:" << exeption_num << ". You shouldn't tap space" << endl;
        }
        else
        {
            cout << "Error code:" << exeption_num << ". You shouldn't write numbers" << endl;
        }
            
        getLetter(); //try one more time
        
    }
    catch (const string exeption)
    {
        cout << "Error code: 5. " << exeption << endl;
        getLetter();
    }
    catch (...) //catch-all
    {
        cout << "There is an exeption with an undetermined type\n";
        getLetter();
    }


    inpt.clear();

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

