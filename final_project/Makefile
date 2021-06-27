CC = gcc
CFLAGS = -std=gnu99 -Wall -O
LDLIBS = -lm


out.txt: main
	./main > out.txt

main.o: main.c funcs.c
	$(CC) $(CFLAGS) -c main.c -o main.o
main: main.o funcs.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)


.PHONEY:
clean:
	$(RM) *.o out* main
