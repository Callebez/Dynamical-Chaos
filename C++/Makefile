
all:
	g++ -Wall  ./source/biffurcation.cpp ./source/main.cpp ./source/penduli.cpp ./source/rungekutta4thSquare.cpp ./source/LyapExp.cpp -o Lytest.o -lm

run:
<<<<<<< HEAD
	g++  -Wall ./source/biffurcation.cpp ./source/main.cpp ./source/penduli.cpp ./source/plotting.cpp ./source/rungekutta4thSquare.cpp -o main.o -lm
	./main.o
=======
	g++  -Wall ./source/biffurcation.cpp ./source/lyapunov.cpp  ./source/LinearAlgebra.cpp ./source/main.cpp ./source/penduli.cpp ./source/lorenz.cpp ./source/rungekutta4thSquare.cpp -o Lytest.o -lm
	./Lytest.o
>>>>>>> Lyapunov-exponent-review
clean:
	find . -type f -executable -exec rm '{}' \;

#   mingw32-make -f Makefile #Para funcionar no Windows, Linux será somente escrever no terminal "make", ou "make run", acho