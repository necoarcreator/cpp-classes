#include "CBinaryTree.h"
#include <iostream>
using namespace std;

template <class T>
CBinaryTree<T>::CBinaryTree() : root{ nullptr } {};

template <class T>
CBinaryTree<T>::CBinaryTree(T num)
{
	Node<T>* perehKnot = new Node<T>(num); //����� ������� �����, �������� ��������� ��� uniform-������������� ���� Node* perehKnot {num}. �������� ������� �� ������� new,��
	root = perehKnot;
}

template <class T>
void CBinaryTree<T>::addKnot(T num)
{
	Node<T>** plceKnot = &root;
	while (*plceKnot != nullptr) {

		Node<T>& tempKnot = **plceKnot;
		if (num < tempKnot.val) {
			plceKnot = &tempKnot.left;
		}
		else if (num > tempKnot.val) {
			plceKnot = &tempKnot.right;
		}
		else {
			cout << "This node already exists!\n";
			return;
		}
	}
	*plceKnot = new Node<T>(num); //����� �������� nullptr ��� ����� ����������

	if ((root->left != nullptr) || (root->right != nullptr))
	{
		cout << "A new node was successefully added!\n";

	}
}

template <class T>
T CBinaryTree<T>::delKnot(T num)
{
	Node<T>** plceKnot = &root;
	while (*plceKnot != nullptr) {

		Node<T>& tempKnot = **plceKnot;
		if (num < tempKnot.val) {
			plceKnot = &tempKnot.left;
		}
		else if (num > tempKnot.val) {
			plceKnot = &tempKnot.right;
		}
		else {

			cout << "Error!\n";
			return 0;
		}
		int flag;
		if (((**plceKnot).val == num) && ((**plceKnot).left == nullptr) && ((**plceKnot).right == nullptr))
		{ //��� ��������
			
			if ((**tempKnot.left) != nullptr)
			{
				flag = 0;
			}
			else if ((**tempKnot.right) != nullptr)
			{
				flag = 1;
			}

				if ((flag == 0) and ((**plceKnot).val == (*tempKnot.left).val)) //������� ����� �����
				{
					T n = (*tempKnot.left).val;
					*plceKnot = nullptr;
					tempKnot.left = nullptr;
					return n;
				}
		
			else if (((**tempKnot.right) != nullptr) and (flag == 1))
			{ 
				if ((**plceKnot).val == (*tempKnot.right).val)//������� ����� ������
				{
				T m = (*tempKnot.right).val;
				*plceKnot = nullptr;
				tempKnot.right = nullptr;
				return m;
				}
			}
			else
			{
				cout << "Error!\n";
			}
		}
		else if ((((**plceKnot).val == num)) && (((**plceKnot).left != nullptr) || ((**plceKnot).right != nullptr)))
		{ //���� ���� ���� �� ���� �������

			if ((**plceKnot).left != nullptr) //���� ���� ����� ������� � ���������� ���� plceKnot
			{
				if ((**plceKnot).val == (*tempKnot.left).val) //������� ����� ����� � tempKnot
				{
					T n = (*tempKnot.left).val;
					tempKnot.left = (**plceKnot).left;
					(**plceKnot).left = nullptr;
					return n;
				}
				else
				{ //������� ����� ������ � tempKnot
					T m = (*tempKnot.right).val;
					tempKnot.right = (**plceKnot).left;
					(**plceKnot).left = nullptr;
					return m;
				}
			}
			else //���� ���� ������ ������� � ���������� ���� plceKnot
			{
				if ((**plceKnot).val == (*tempKnot.left).val) //������� ����� ����� � tempKnot
				{
					T n = (*tempKnot.left).val;
					tempKnot.left = (**plceKnot).right;
					(**plceKnot).right = nullptr;
					return n;
				}
				else //������� ����� ������ � tempKnot
				{
					T m = (*tempKnot.right).val;
					tempKnot.right = (**plceKnot).right;
					(**plceKnot).right = nullptr;
					return m;
				}
			}
		}

		else if (((**plceKnot).val == num) and (((**plceKnot).left != NULL) and ((**plceKnot).right != NULL))) //������� ����� ���� � ����� ������������
		{
			T n = (**plceKnot).val;
			Node<T>& childKnot = *(**plceKnot).right; //�������������
			Node<T>& fatherKnot = **plceKnot; //������ �������������. �������� � ������� ��������� ���������� ����

			while ((childKnot).left != nullptr) //���� �� ����� �� ������������ (�������)
			{
				fatherKnot = childKnot;
				childKnot = *childKnot.left;

			}
			
			(**plceKnot).val = (childKnot).val;
			(fatherKnot).left = nullptr;
			return n;
		}
	}
	*plceKnot = new Node<T>(num); //����� �������� nullptr ��� ����� ����������

}


template <class T>
Node<T>* CBinaryTree<T>::getroot()
{
	return root;
}

template <class T>
void CBinaryTree<T>::print(Node<T>* root) {
	if (root != nullptr) {

		print(root->left);

		cout << root->val << " -> ";

		print(root->right);
	}
}

template <class T>
T CBinaryTree<T>::SrchPlc(T num)
{
	T counter = 0;
	Node<T>** plceKnot = &root;
	while (*plceKnot != nullptr) {
		counter++;
		Node<T>& tempKnot = **plceKnot;
		if (num < tempKnot.val) {
			plceKnot = &tempKnot.left;
		}
		else if (num > tempKnot.val) {
			plceKnot = &tempKnot.right;
		}
		else {

			cout << "Error!\n";
			return 0;
		}

		if ((**plceKnot).val == num)
		{
			if ((**plceKnot).val == (*tempKnot.left).val)
			{
				T n = (*tempKnot.left).val;
				cout << "It is " << counter << "th left node.\n";
				return n;
			}
			else
			{
				T m = (*tempKnot.right).val;
				cout << num << " is " << counter << "th right node.\n";
				return m;
			}
		}
	}
	*plceKnot = new Node<T>(num); //����� �������� nullptr ��� ����� ����������

}

