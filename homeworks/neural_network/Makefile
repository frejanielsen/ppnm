CC = gcc
DIR = ../minimum
CFLAGS += -I$(DIR)
CFLAGS += -Wall -std=c11 $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)


.PNONEY: default
#default: output_ann.txt plot.png

default: output.txt
	cat output.txt

output.txt: main
	./main > output.txt

main: main.o ann_fs.o

main.o: main.c ann_fs.h

ann_fs.o: ann_fs.c
	$(CC) $(CFLAGS) $(CPPFLAGD) -o $@ -c $<


#plot.png: output_ann.txt Makefile
#	echo '\
#	set terminal png;\
#	set output "$@";\
#	set key top left;\
#	set tics out;\
#	set title "Neural networks interpolarization";\
#	plot "$<" index 0 using 1:2 with lines title "interpolated data"\
#	,"$<" index 1 using 1:2 with points title "data"\
#	' | gnuplot

.PHONEY: clean
clean:
	$(RM) *.o main out.*

.PHONEY: test
test:
	echo $(LDLIBS)
	echo $(CC)
	echo $(RM)
