
# Compiler settings - Can be customized.
CC = g++

all: clean run.exe

run.exe: Queue.o Stack.o lifofifo_template.o
	$(CC) lifofifo_template.o Queue.o Stack.o -o run.exe

lifofifo_template.o: lifofifo_template.cpp
	$(CC) -c lifofifo_template.cpp

Stack.o: Stack.cpp
	$(CC) -c Stack.cpp

Queue.o: Queue.cpp
	$(CC) -c Queue.cpp

clean:
	rm -rf ./target/*.out