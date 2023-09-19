

#include <iostream>
#include <vector>
#include <memory>
#include <string>

int main()
 {
    std::string s, words;
    int i, j;

    std::cin >> s; // считываем строку

    std::vector <char> word{s.begin(), s.end()}; //записываем посимвольно вспомогательный вектор-слово

    j = word.size(); //устнавливаем число попыток

    
    
    std::cout << "The game has started. Guess a letter!" << std::endl;
    std::string ch(j, '-'); //это строка, которую будем выводить каждый раз на экран
    
    std::cout << ch << std::endl;

    char buk;
    int k;
    bool FL = false; //флаг, поможет запускать различные блоки кода

    for (i = 0; i < j; i++) 
    {
        
        std::cin >> buk;
        for (k = 0; k < j; k++) //блок проверки на наличие буквы в слове
        {
            if (word[k] == buk) 
            {
                ch[k] = buk;
                FL = true;
                
                
            }
        }
        if (FL == true) {
            i -= 1; //если хотя бы одна буква угадана, то жизнь "сохраняется"
            std::cout << "Correct guess!" << std::endl;
        }
        else {
            std::cout << "Incorrect..." << std::endl;
        }

        std::cout << ch << std::endl; //готовимся к следующей попытке

        FL = false; 


        if (ch == s)
        {
            std::cout << "You won! The word was:" << std::endl;
            std::cout << ch << std::endl;

            FL = true; 
            break;
        }
    }
    if (FL == false) {
        std::cout << "You've lost! The word was:" << std::endl;
        std::cout << s << std::endl;
    }
    
    return 0;
}

/* {
    std::unique_ptr<int> sm = std::make_unique<int>(17);
    std::cout << *sm << std::endl;

*/















