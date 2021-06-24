CFLAGS = -Wall -Werror  -O -std=gnu11 #-DNDEBUG
LDLIBS = -lm


.PHONEY: default
default: out.main.txt
	cat  "$<"

out.main.txt: main
	./$< > $@


main: main.o mc_fs.o

PHONEY: clean

clean:
	$(RM) out* main *.o
