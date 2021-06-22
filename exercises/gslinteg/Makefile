CFLGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

.PHONEY: all
all: out.value_exercise_A.txt

out.value_exercise_A.txt: main Makefile
	./$<

.PHONEY: clean
clean:
	$(RM) main out* *.png
