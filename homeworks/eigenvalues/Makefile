CFLAGS = -Wall  -O -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


.PHONEY: default


default:  out_ev.txt States.png
	cat out_ev.txt

out_ev.txt: main
	./main > out_ev.txt

main.o: main.c evfs.c

main: main.o evfs.o

States.png: out_states.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "State";\
	set key bottom right;\
	Pi=3.1415926535897932348;\
	a1=1;\
	a2=1;\
	a3=1;\
	S1(x)=a1*sin(Pi*x);\
	S2(x)=a2*sin(Pi*2*x);\
	S3(x)=a3*sin(Pi*3*x);\
	fit S1(x) "$<" index 0 via a1;\
	fit S2(x) "$<" index 1 via a2;\
	fit S3(x) "$<" index 2 via a3;\
	plot "$<" index 0 using 1:2 w lp title "numerical state 1"\
	, "$<" index 1 using 1:2 w lp title "numerical state 2"\
	, "$<" index 2 using 1:2 w lp title "numerical state 3"\
	, S1(x) w l title "exact state 1 norm fitted"\
	, S2(x) w l title "exact state 2 norm fitted"\
	, S3(x) w l title "exact state 3 norm fitted";\
	' | gnuplot

.PHONEY: clean

clean:
	$(RM) out* main *.o *.png *.log

