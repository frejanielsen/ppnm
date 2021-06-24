CFLAGS = -Wall -Werror -O -std=gnu11
LDLIBS = -lm
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: out.main.txt
	cat  out.main.txt

out.main.txt: main
	./$< > $@

main: main.o funcs.o


.PHONEY: clean

clean:
	$(RM) out* main *.o *.png *.log
