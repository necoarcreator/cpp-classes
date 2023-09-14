
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

int first(void);
int sec(void);
int thrd(void);
int frth(void);
int ffth(void);
int sxth(void);

class Point
{
public:
    Point()
    {
    }
    Point(int, int, int, int, int)
    {
    }
};

int main(void)
{ 
    first();
    
    sec();

    thrd();

    printf("почему-то сразу 4 функции не запусккаются");

    sxth();

    return 0;
}


int first(void)
{
    int m[5];
    int i;
    for (i = 0; i < 5; i++)
    {
        m[i] = i;
        printf("m[%1d] = %2d\n", i, m[i]);
    }
    return 0;
}

int sec(void)
{
    int* p;
    int i, n;
    p = (int*)malloc(sizeof(int));
    printf("Enter the number of elements: ");
    // 'scanf' is deprecated: This function or variable may be unsafe. Consider using scanf_s instead. //
    scanf_s("%d", &n);
    if (p == NULL)
    {
        printf("process executed with flag 1");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < n; ++i)
    {
        p[i] = i;
        printf("p[%1d] = %2d\n", i, p[i]);
    }
    free(p);
    return 0;
}
int thrd(void)
{
    int m[5]{ 0, 1, 2, 3, 4 };
    int i;
    for (i = 0; i < 5; i++) 
    {
        printf("m[%1d] = %2d\n",i, m[i]);
    }
}

//int frth(void)
//{   //тут не особо что-то получилось. было переполнение данных, хотя задать что-то классом звучало заманчиво...
    //Point m[5] = {
    //    Point(0, 1, 2, 3, 4)
    //};

    //int i;
    //for (i = 0; i < 5; i++)
    //{
    //    printf("m[%1d] = %2d\n", i, m[i]);
    //}

//}

int ffth(void)
{   
    // не совсем понимаю, как из массива указателей вытащитть значения
    int u = 0, v = 1, g = 2, h = 3, c = 4;

    int* const rr[5] = { &u, &v, &g, &h, &c };
    int i;

    for (i = 0; i < 5; i++)
    {
        std::cout << rr[i] << std::endl;
    }



}

int sxth(void)
{
    int* p, * q, i;
    p = (int*)calloc(3, sizeof(int));
    if (p == NULL)
    {
        std::cout << "Не удалось выделить память" << std::endl;
        exit(EXIT_FAILURE);
    }

    q = (int*)realloc(p, 5 * sizeof(int));
    if (q == NULL)
    {
        std::cout << "Не удалось увеличить память" << std::endl;
    }
    else {
        p = q;
        for (i = 0; i < 5; i++) {
            p[i] = i;
            printf("p[%1d] = %2d\n", i, p[i]);
        }

    }
    free(p);
}
    
//a[n]//