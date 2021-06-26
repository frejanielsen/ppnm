CFLAGS = -Wall -Werror  -O -std=gnu11 #-DNDEBUG
LDLIBS = -lm


.PHONEY: default
default: out.main.txt out_err.png
	cat  "$<"

out.main.txt: main
	./$< > $@


main: main.o mc_fs.o


ERRPLOT =  out_err.png
ERRDATA =  out_err.txt
ERRGPI  =  err.gpi
$(ERRPLOT): $(ERRGPI) $(ERRDATA)
	pyxplot $<




PHONEY: clean

clean:
	$(RM) out* main *.o
