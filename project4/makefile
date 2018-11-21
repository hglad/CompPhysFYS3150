CC = gcc
CXX = mpic++
CPPflags = -g -Wall -std=c++11 -o3

INCL :=
LIB := -DARMA_DONT_USE_WRAPPER -lblas -llapack


main.x: main_ising_parallel.o ising.o
	${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}
