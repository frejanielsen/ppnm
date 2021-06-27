CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: out.txt
	cat out.txt

out.txt: main
	./main > out.txt

main.o: main.c funcs.c
	$(CC) $(CFLAGS) -c main.c -o main.o
main: main.o funcs.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)


.PHONEY: clean
clean:
	$(RM) *.o out* main
