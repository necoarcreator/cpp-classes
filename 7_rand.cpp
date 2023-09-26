#include <iostream>
#include <io.h>
#include <fcntl.h>
#include <string>
#include <fstream>
using namespace std;

int main()
{
    _setmode(_fileno(stdout), _O_U16TEXT);
    _setmode(_fileno(stdin), _O_U16TEXT);
    _setmode(_fileno(stderr), _O_U16TEXT);
    wstring userName, birthDate;
    int i, sum = 0;
    wcout << L"Введите своё имя:" << endl;

    getline(wcin, userName);

    wcout << L"Введите свою дату рождения в формате ЧЧ.ММ.ГГГГ:" << endl;

    getline(wcin, birthDate);

    for (i = 0; i < birthDate.length(); i++)
    {
        if ((i != 2) && (i != 5))
        {
            sum += static_cast<int>(birthDate[i]) - static_cast<int>('0');

        }


    }

    const int monthPlaceStart = 3;

    const int monthPlaceLength = 2;

    wstring month = birthDate.substr(monthPlaceStart, monthPlaceLength);

    int j, sumMonth = 0;
    for (j = 0; j < month.length(); j++)
    {
        sumMonth += static_cast<int>(month[j]) - static_cast<int>('0');

    }

    int randPool[3]{0, 11111111, 22222222}, randNum;

    srand(time(NULL));

    randNum = randPool[0 + rand() % (2 - 0 + 1)];

    randNum += sum ^ (sumMonth);

    int numMusic = 0b1111, numFilm = 0b1111, numActor = 0b1111, numActress = 0b1111, numCelebrity = 0b1111;

    numMusic = numMusic & randNum;
    randNum = randNum >> 4;

    numFilm = numFilm & randNum;
    randNum = randNum >> 4;

    numActor = numActor & randNum;
    randNum = randNum >> 4;

    numActress = numActress & randNum;
    randNum = randNum >> 4;

    numCelebrity = numCelebrity & randNum;
    randNum = randNum >> 4;


    wstring music[14]{ L"Master of puppets Metallica", L"Брошу Pyrokinesis",
    L"Homebody PH1", L"Червяк Даргомыжский", L"Russian Ebunny",
    L"Enemy Imagine dragons", L" Зизазай Огги и тараканы",
    L"Valentine Justice", L"Swimming pools Kendrick Lamar",
    L"Night Call Kavinsky", L"Yesterday The Beatles",
    L"Так закалялась сталь Гражданская оборона",
    L"За деньги да Инстасамка", L"Bad habits Ed sheeran" };
    
    wstring movies[14]{ L"Короткое замыкание", L"Трансформеры", L"Полночь",
    L"По соображениям совести", L"Джентльмены", L"Ван Хельсинк",
    L"Операция Ы", L"Драйв", L"Большая игра", L"Игра в имитацию",
    L"Всё везде и сразу", L"Форрест Гамп", L"Терминал",
    L"Не смотрите наверх" };

    wstring actors[14]{ L"Райан Гослинг", L"Кристиан Бейл", L"Джеки Чан",
    L"Эндрю Гарфилд", L"Леонардо ДиКаприо", L"Бенедикт Кембербетч",
    L"Мэтью Макконахи", L"Дэвид Линч", L"Дэниэлл Рэдклифф",
    L"Джейк Джиллехолл", L"Александр Петров", L"Рональд Рейган",
    L"Том Круз", L"Брэд Питт" };

    wstring actriss[14]{ L"Наталия Крачковская", L"Меган Фокс", L"Джениффер Лопез",
    L"Эмма Стоун", L"Эмма Уотсон", L"Кира Найтли", L"Милла Йовович",
    L"Марго Робби", L"Мерил Стрип", L"Скарлет Йоханссон", L"Александра Бортич",
    L"Анджелина Джоли", L"Шарлиз Терон", L"Ингеборга Дапкунайте" };

    wstring celeb[14]{ L"Дейв Майнстейн", L"Иван Зола", L"Сергей Мавроди",
    L"Екатерина Гордеева", L"Роберт Опенгеймер", L"Лионель Месси",
    L"Папич", L"Хидэо Кодзима", L"Ляйсан Утяшева", L"Павел Воля", L"Криштиану Роналду",
    L"Борис Бурда", L"Ким Кардашьян", L"Клава Кока" };
        
    wcout << L"Добрый день, " << userName << L", ваш рандом на сегодня!" << endl;

    wcout << L"Вашa песня на сегодня: " << music[numMusic] << endl;


    wcout << L"Ваш фильм: " << movies[numFilm] << endl;



    wcout << L"Ваш актер: " << actors[numActor] << endl;



    wcout << L"Ваша актриса: " << actriss[numActress] << endl;



    wcout << L"Ваша селебрити: " << celeb[numCelebrity] << endl;
}
