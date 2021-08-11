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

