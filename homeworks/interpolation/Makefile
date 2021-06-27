CC = gcc
CFLAGS = -O -Werror -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: out.lin.txt out.quad.txt out.cube.txt out.gsl.lin.txt interp.png integ.png integ.quad.png


out.lin.txt: lin
	./$< > $@

out.quad.txt: quad
	./$< > $@

out.cube.txt: quad
	./$< > $@

out.gsl.lin.txt: gsl_linterp
	./$< > $@

interp.png: out.lin.txt out.quad.txt out.cube.txt out.xydata.txt out.gsl.lin.txt out.gsl.pol.txt Makefile 
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set key top;\
	plot "$<" using 1:2 w l title "lin interp"\
	,    "out.quad.txt" w l title "quad interp"\
	,    "out.cube.txt" w l title "cubic interp"\
	,    "out.xydata.txt" using 1:2 w p title "data"\
	,    "out.gsl.lin.txt" using 1:2 w l title "gsl linear interp"\
	,    "out.gsl.pol.txt" using 1:2 w l title "gsl polynomial interp"\
	' | gnuplot



integ.png: out.xyinteg.txt out.integ.exact.txt out.xyinteg.half.txt out.gsl.lin.integ.txt Makefile 
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set key bottom left;\
	plot "$<" using 1:2 w l title "lin integ"\
	,    "out.integ.exact.txt" using 1:2 w p title "analytical integral"\
	,    "out.xyinteg.half.txt" using 1:2 w l title "lin integ using half the points"\
	,    "out.gsl.lin.integ.txt" using 1:2 w l title "gsl linear integration"\
	' | gnuplot



integ.quad.png: out.quad.integ.txt out.integ.exact.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set key bottom left;\
	plot "$<" using 1:2 w l title "quad integ"\
	,    "out.integ.exact.txt" using 1:2 w p title "analytical integral"\
	' | gnuplot




.PHONEY: clean
clean:
	$(RM) *.o lin out*

.PHONEY: test
test:
	echo $(LDLIBS)
	echo $(CC)
	echo $(RM)
