#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <chrono>
#include <locale>
#include <codecvt>
#include <sstream>

#include <iomanip>
#include <windows.h>

using namespace std;


class dictionary
{
    const locale utf8_locale = locale(locale(), new codecvt_utf8<wchar_t>());
    int random = 0, length = 0;
    vector<wstring> dict;
    wstring wrword = L"";
    wstringstream rdword;


public:

    wstring answer = L"n";

    void setDict(int num = 1)
    {
        wifstream buf(L"dict1.txt");
        buf.imbue(utf8_locale);
        if (!buf)
        {
            cerr <<"О нет, dict1.txt нельзя открыть для чтения!" << endl;
            exit(1);
        }


        while (buf)
        {
            rdword << buf.rdbuf();
            length += 1;
            dict.push_back(rdword.str());
        }
        buf.close();
        rdword.clear();
    }

    wstring getWord()
    {
        setDict(1);
        srand(time(NULL));
        this_thread::sleep_for(chrono::seconds(1));
        random = rand() % length;
        return dict[random];
    }
    int addWord()
    {

        wcout << "Хотите добавить новое слово? Нажмите [Y/n]" << endl;
        wcin >> answer;
        if (answer == L"Y")
        {
            wofstream bufIn(L"dict1.txt", ios::app);
            bufIn.imbue(utf8_locale);
            wcout << "Ну так напишите его скорее: ";
            wcin >> wrword;
            bufIn << endl;
            bufIn << wrword;
            bufIn.close();
        }
        else
        {
            return 0;
        }
        addWord();
    }

    
};

int main()
    {
    SetConsoleOutputCP(1251);
    SetConsoleCP(1251);
        dictionary dict;
       dict.addWord();
    }