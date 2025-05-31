#!/bin/bash
cp exampleRequest 
./execute calc "16265 + 889219.029"

echo "+ operation OK"

./execute calc "16265 + 889219.029 - 26267.27"

echo "- operation OK"

./execute calc "16265 + 889219.029 - (26267.27 * 21 )"

echo "* operation OK"


./execute calc "16265 + 889219.029 - (26267.27 * 21) / 2"

echo "/ operation OK"

./execute calc "16265 + 889219.029 - (26267.27 * 21 - 145 % 17) / 2"

echo "% operation OK"

./execute calc exampleSingleRequest

echo "One-string file check OK"

./execute calc exampleRequest

echo "Multistring file check OK"

./execute calc -l "16265 + 889219.029 - (26267.27 * 21 - 145 % 17) / 2"

echo "shell output OK"

./execute calc -l log testOutput "16265 + 889219.029 - (26267.27 * 21 - 145 % 17) / 2"

echo "log output OK"

./execute calc -l log testOutput exampleSingleRequest

echo "file read and write OK"

echo "All tests are successful"

