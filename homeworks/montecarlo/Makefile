CFLAGS = -Wall -Werror  -O -std=gnu11 #-DNDEBUG
LDLIBS = -lm


.PHONEY: default
default: out.main.txt err.png
	cat  "$<"

out.main.txt: main
	./$< > $@


main: main.o mc_fs.o

err.png: out_err.txt Makefile
	echo '\
	set terminal png;\
	set output "err.png";\
	set xlabel "Number of points";\
	set ylabel "Error";\
	set key top;\
	set format x "%.0e";\
	set format y "%.0e";\
#	a=1; b=-0.5; f(x)=a*x**b;\
#	c=1; d=-1; g(x)=c*x**d;\
#	fit f(x) "$<" using 1:2 via a,b;\
#	fit g(x) "$<" using 1:3 via c,d;\
	plot "$<" using 1:2 w p title "multi monte carlo"\
	, "$<" using 1:3 w p title "plain monte carlo"\
#	, f(x) w l title "fit1"\
#	, g(x) w l title "fit2";\
	' | gnuplot




PHONEY: clean

clean:
	$(RM) out* main *.o
