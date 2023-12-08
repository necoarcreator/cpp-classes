#include "Node.h"
using namespace std;

template <class T>

Node<T>::Node() : val{ 0 }, left{ nullptr }, right{ nullptr } {};

template <class T>

Node<T>::Node(T num) : val{ num }, left{ nullptr }, right{ nullptr } {};

