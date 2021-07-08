
all:
	g++ -Wall  ./source/biffurcation.cpp ./source/main.cpp ./source/penduli.cpp ./source/rungekutta4thSquare.cpp -o   main.o -lm

run:
	g++  -Wall ./source/biffurcation.cpp ./source/main.cpp ./source/penduli.cpp ./source/rungekutta4thSquare.cpp -o main.o -lm
	./main
clean:
	find . -type f -executable -exec rm '{}' \;

#   mingw32-make -f Makefile #Para funcionar no Windows, Linux será somente escrever no terminal "make", ou "make run", acho