#!/bin/bash

# Удаление файлов в поддиректориях csv
for dir in ./out/csv/{1,2,3,test}/{1..8}; do
    if [ -d "$dir" ]; then
        rm -f "$dir"/*
        #echo "Удалены файлы в $dir"
    else
        #echo "Директория $dir не существует, пропускаем"
    fi
done

# Удаление файлов в директории pics
if [ -d "./out/pics" ]; then
    rm -f ./out/pics/*
    #echo "Удалены файлы в ./out/pics"
else
    #echo "Директория ./out/pics не существует, пропускаем"
fi

#echo "Очистка завершена"