CFLAGS = -Wall -O -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
#LDLIBS = -lm
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default

default: out.txt
	cat out.txt


out.txt: main
	./main > out.txt

main.o: main.c min_fs.h

main: main.o min_fs.o

.PHONEY: clean

clean:
	$(RM) out* main *.o
