CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default

default: out.txt ls_plot.png
	cat out.txt


ls_plot.png:  out_fit.txt out_data.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "Time (days)";\
	set ylabel "Counts";\
	set key top;\
	plot "$<" index 0 using 1:2 w l title "best fit"\
	, "$<" index 1 using 1:2 w l title "fit + Delta c"\
	, "$<" index 2 using 1:2 w l title "fit - Delta c"\
	, "out_data.txt" using 1:2:3 with errorbars title "Data";\
	' | gnuplot

out.txt: main
	./main > out.txt

main.o: main.c ls_fs.c

main: main.o ls_fs.o

.PHONEY: clean
clean:
	$(RM) out.* *.o main *.png
