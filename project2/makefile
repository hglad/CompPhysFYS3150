CPPflags = g++ -O3

LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack

2d: main_d.o project2.o
	${CPPflags} -o project2d.x main_d.cpp project2.cpp ${LIB}
