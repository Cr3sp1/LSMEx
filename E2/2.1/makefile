CC = g++
CFLAGS = -Wall -O3 --std=c++11
LDFLAGS = -larmadillo 

main.exe : main.o random.o LSMfunctions.h Functions.h Integral.h
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS) $(LDFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
run: 
	./main.exe
clean :
	rm *.o main.exe seed.out
