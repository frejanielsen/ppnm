CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)



.PHONEY: default

default: out_lin.txt
	cat out_lin.txt

out_lin.txt: main
	./main > out_lin.txt

main.o: main.c

main: main.o

.PHONEY: clean
clean:
	$(RM) out.* *.o main *.png
