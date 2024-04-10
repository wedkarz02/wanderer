CPPFLAGS=-std=c++11 -fPIE -O3 -g

main: matrix.o main.o sparse.o
	g++ -o runme matrix.o sparse.o main.o $(CPPFLAGS)
matrix.o: matrix.cpp matrix.h
	g++ -c matrix.cpp $(CPPFLAGS)
sparse.o: sparse.cpp sparse.h
	g++ -c sparse.cpp $(CPPFLAGS)
main.o: main.cpp
	g++ -c main.cpp  $(CPPFLAGS)

clean:
	rm *.o runme
