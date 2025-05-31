
all: run

run: main.o krootoy.o godunov2.o tests.o
	g++ main.o krootoy.o godunov2.o tests.o -o run
	
main.o: main.cpp
	g++ -c main.cpp
	
krootoy.o: krootoy.cpp tests.o
	g++ -c krootoy.cpp
	
godunov2.o: godunov2.cpp krootoy.o tests.o
	g++ -c godunov2.cpp

tests.o: tests.cpp
	g++ -c tests.cpp

clean:
	rm -rf *.o run