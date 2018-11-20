MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
MPI_LINK_FLAGS = $(shell mpicc --showme:link)

#CPPflags = g++ -O3

LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack

main: main_ising_parallel.o ising.o
	${CC} $(MPI_COMPILE_FLAGS) main_ising_parallel.cpp ising.cpp $(MPI_LINK_FLAGS) -o main.x ${LIB
