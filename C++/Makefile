CC := g++
CFLAGS := -g -Wall -lm -pthread

all:
	$(CC) $(CFLAGS)  ./source/biffurcation.cpp ./source/lyapunov.cpp ./source/plotting.cpp ./source/printing.cpp ./source/LinearAlgebra.cpp ./source/main.cpp ./source/penduli.cpp ./source/lorenz.cpp ./source/rungekutta4thSquare.cpp ./source/discretelyap.cpp -o Lytest.o 

run: 
	$(CC) $(CFLAGS)  ./source/biffurcation.cpp ./source/lyapunov.cpp ./source/plotting.cpp ./source/printing.cpp ./source/LinearAlgebra.cpp ./source/main.cpp ./source/penduli.cpp ./source/lorenz.cpp ./source/rungekutta4thSquare.cpp ./source/discretelyap.cpp -o Lytest.o 
	 ./Lytest.o
clean:
	find . -type f -executable -exec rm '{}' \;

#   mingw32-make -f Makefile #Para funcionar no Windows, Linux será somente escrever no terminal "make", ou "make run", acho