#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <chrono>

using namespace std;

class dictionary
{
    int random = 0, length = 0;
    vector<string> dict;
    string word = "";

public:

    string answer = "n";

    void setDict(int num = 1)
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

    string getWord()
    {
        setDict(1);
        srand(time(NULL));
        this_thread::sleep_for(chrono::seconds(1));
        random = rand() % length;
        return dict[random];
    }
    int addWord()
    {

        cout << "Do you want to add the new word? Tap [Y/n]" << endl;
        cin >> answer;
        if (answer == "Y")
        {
            cout << "You shoould write it in lower register." << endl;
            fstream bufIn("dict1.txt", ios::app);
            cout << "Insert new word: ";
            cin >> word;
            bufIn << "\n"; //почему-то два отступа
            bufIn << word;
            bufIn.close(); //не знаю, как пофиксить следующее: добавление слова ведёт к вылету проги
        }
        else
        {
            return 0;
        }
        addWord();
    }
    ~dictionary() {}
};

/*
int main()
{
    Dictionary dict;
    cout << dict.getWord() << endl;
    dict.addWord();
}
*/
