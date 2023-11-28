#ifndef _IRF_H
#define _IRF_H


class Interface
{
public:
    virtual void push(int num) = 0;
    virtual int pop() = 0;
    virtual bool isEmpty() = 0;
};

#endif