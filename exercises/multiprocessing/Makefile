CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += -pthread
CFLAGS += -fopenmp
LDLIBS = -lm -lgomp -lpthread

OBJECTS = main.o random_number.o pi_approx.o
HDRS = random_number.h pi_approx.h
EXECS = main

default: $(EXECS)
	./$<

$(EXECS): $(OBJECTS)

$(OBJECTS): $(HDRS)

clean:
	$(RM) $(EXECS) $(OBJECTS)
