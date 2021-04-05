CC =  gcc
CFLAGS = -O -fwrapv -Wall
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: epsilon
	./epsilon > out.txt

epsilon: epsilon.o
	$(CC) -o epsilon epsilon.o $(LDLIBS)

epsilon.o: epsilon.c
	$(CC) $(CFLAGS) -c epsilon.c

clean:
	$(RM) epsilon.o epsilon out.txt

#.PHONEY: test
#test: 			#test target used for debugging
#	echo $(LDLIBS)
#	echo $(CC)
#	echo $(RM)
