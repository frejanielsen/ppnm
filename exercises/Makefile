CFLAGS = -O -W -std=gnu99
CFLAGS += $(gsl-config --cflags)
LDLIBS += $(gsl-config --libs)

all: output.txt ; cat $<
output.txt: main ; ./$< > $@
clean: ; $(RM) main output
