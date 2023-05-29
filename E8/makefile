CC = g++
CFLAGS = -Wall -O3 --std=c++11 
LDFLAGS = -larmadillo

main.exe : main.o random.o LSMfunctions.h lez8.h
	$(CC) random.o main.o -o main.exe $(CFLAGS) $(LDFLAGS)
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS) $(LDFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS) $(LDFLAGS)
run: 
	./main.exe
clean :
	rm *.o main.exe seed.out
