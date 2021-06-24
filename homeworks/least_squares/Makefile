CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: clean all

out.txt: main
	./main > out.txt

main.o: main.c ls_fs.c

main: main.o ls_fs.o

clean:
	$(RM) out.* *.o main *.png
