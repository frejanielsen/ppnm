CC = gcc
CFLAGS = -O  -std=gnu11
LDLIBS = -lm


exp.pdf: exp.tex exp_plot.pdf
	pdflatex $<
	pdflatex $<

out.txt:main
	./$< > $@

main:main.o

exp_plot.pdf: exp_plot.gpi out.txt
	pyxplot $<

clean:
	$(RM) exp.pdf *.log *.aux out* main *.o fig*.pdf fig*.tex
