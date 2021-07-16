
all:
	gcc -Wall   ./source/2sp-osc.c ./source/cgic.c ./source/control.c ./source/dverk.c ./source/dynamo.c ./source/main.c ./source/rinaldi.c ./source/rossler.c  -o DverkTest.o -lm

run:
	gcc  -Wall  ./source/2sp-osc.c ./source/cgic.c ./source/control.c ./source/dverk.c ./source/dynamo.c ./source/main.c ./source/rinaldi.c ./source/rossler.c  -o DverkTest.o -lm
	./DverkTest.o
#	g++  -Wall ./source/biffurcation.cpp ./source/main.cpp ./source/penduli.cpp ./source/rungekutta4thSquare.cpp ./source/LyapExp.cpp -o Lytest.o -lm
#	./Lytest.o
clean:
	find . -type f -executable -exec rm '{}' \;

#   mingw32-make -f Makefile #Para funcionar no Windows, Linux será somente escrever no terminal "make", ou "make run", acho