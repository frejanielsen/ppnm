CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: out_root.txt
	cat out_root.txt

out_root.txt: main
	./main > out_root.txt

main.o: main.c rootfs.c linfs.c odefs.c

main: main.o rootfs.o linfs.o odefs.o

.PHONEY:clean

clean:
	$(RM) out.* *.o main *.png
