CC = g++
CFLAGS = -Wall

All: execute

execute: main.o Helper.o StandartCalc.o FileWriteCalc.o FileReadCalc.o FileReadWriteCalc.o
	$(CC) $(CFLAGS) main.o Helper.o StandartCalc.o FileWriteCalc.o FileReadCalc.o FileReadWriteCalc.o -o execute

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp -o main.o

Helper.o: Helper.cpp
	$(CC) $(CFLAGS) -c Helper.cpp -o Helper.o

StandartCalc.o: StandartCalc.cpp
	$(CC) $(CFLAGS) -c StandartCalc.cpp -o StandartCalc.o

FileWriteCalc.o: FileWriteCalc.cpp
	$(CC) $(CFLAGS) -c FileWriteCalc.cpp -o FileWriteCalc.o

FileReadWriteCalc.o: FileReadWriteCalc.cpp
	$(CC) $(CFLAGS) -c FileReadWriteCalc.cpp -o FileReadWriteCalc.o

FileReadCalc.o: FileReadCalc.cpp
	$(CC) $(CFLAGS) -c FileReadCalc.cpp -o FileReadCalc.o


clean: 
	@echo "Cleaning up..."
	rm execute
	rm -rf *.o execute 
	

rebuild: clean execute

run: execute
	./execute:
