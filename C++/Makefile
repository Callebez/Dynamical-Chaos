target 	:= main 
OBJDIR  := ./build/ 
objs 	:=  ./build/main.o ./build/lyapunov.o ./build/biffurcation.o ./build/plotting.o ./build/printing.o ./build/LinearAlgebra.o ./build/penduli.o ./build/penduli.o ./build/lorenz.o ./build/rungekutta4thSquare.o
CC 		:= g++
CFLAGS 	:= -g -Wall -lm 

all: $(target)

deps := $(patsubst %.o, %.d, $(objs))
# -include $(deps)
DEPFLAGS = -MMD -MF $(@:.o=.d)

main: $(objs)
	$(CC) $(CFLAGS) -o $@ $^

build/%.o: ./source/%.cpp	
	$(CC) $(CFLAGS) -c  $<  -o $@ $(DEPFLAGS)

clean:
	rm -f $(target) $(objs) $(deps)

# all:
# 	g++ -Wall  ./source/biffurcation.cpp ./source/main.cpp ./source/systems.cpp ./source/rungekutta4thSquare.cpp ./source/LyapExp.cpp ./source/FastLyap.cpp -o Lytest.o -lm

# run: 
# 	$(CC) $(CFLAGS)  ./source/biffurcation.cpp ./source/lyapunov.cpp ./source/plotting.cpp ./source/printing.cpp ./source/LinearAlgebra.cpp ./source/main.cpp ./source/penduli.cpp ./source/lorenz.cpp ./source/rungekutta4thSquare.cpp -o Lytest.o 
# 	 ./Lytest.o
# clean:
# 	find . -type f -executable -exec rm '{}' \;

#   mingw32-make -f Makefile #Para funcionar no Windows, Linux será somente escrever no terminal "make", ou "make run", acho