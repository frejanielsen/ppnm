CC =  gcc
CFLAGS = -O -fwrapv -std=gnu11 -Wall
LDFLAGS = -I.
LDLIBS = -lm



default: out.txt
	cat out.txt

out.txt: epsilon  equals
	./$< > $@
	./equals >> $@

.PHONEY: clean
clean:
	$(RM) epsilon.o epsilon out.txt

.PHONEY: test

test:
	echo $(CC)
	echo $(CDLIBS)
	echo $(LDFLAGS)
