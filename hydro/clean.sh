#!/bin/bash

# �������� ������ � �������������� csv
for dir in ./out/csv/{1,2,3,test}/{1..8}; do
    if [ -d "$dir" ]; then
        rm -f "$dir"/*
        #echo "������� ����� � $dir"
    else
        #echo "���������� $dir �� ����������, ����������"
    fi
done

# �������� ������ � ���������� pics
if [ -d "./out/pics" ]; then
    rm -f ./out/pics/*
    #echo "������� ����� � ./out/pics"
else
    #echo "���������� ./out/pics �� ����������, ����������"
fi

#echo "������� ���������"