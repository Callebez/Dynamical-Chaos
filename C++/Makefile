
all:
	g++ -Wall  ./source/biffurcation.cpp ./source/main.cpp ./source/systems.cpp ./source/rungekutta4thSquare.cpp ./source/LyapExp.cpp ./source/FastLyap.cpp -o Lytest.o -lm

run:
	g++  -Wall ./source/biffurcation.cpp ./source/lyapunov.cpp ./source/plotting.cpp ./source/printing.cpp ./source/LinearAlgebra.cpp ./source/main.cpp ./source/penduli.cpp ./source/lorenz.cpp ./source/rungekutta4thSquare.cpp -o Lytest.o -lm
	./Lytest.o
clean:
	find . -type f -executable -exec rm '{}' \;

#   mingw32-make -f Makefile #Para funcionar no Windows, Linux será somente escrever no terminal "make", ou "make run", acho