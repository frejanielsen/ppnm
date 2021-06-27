CC = gcc
CFLAGS =  -O1 -std=gnu11 -Wall
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: result.txt output.txt #plot.png

output.txt: main
	cat result

result.txt: main
	./$< > $@
	cat result.txt

main: main.o ann_fs.o ann_fs.h
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c ann_fs.h
	$(CC) $(CFLAGS) -c $< -o $@

ann_fs.o: ann_fs.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) *.txt *.o main Results.txt
