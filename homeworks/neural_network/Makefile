CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat
PLTUTIL =  pyxplot

EXE   	=  main
OBJ  	=  main.o ann_fs.o input.o
HDRS    =  ann_fs.h input.h
TXT     =  output.txt
DATA    =  data.txt output_ann.txt
GPI     =  plot.gpi
IMGFL   =  plot.png

.PHONEY: default
default: print $(IMGFL)

$(IMGFL): $(GPI) output_ann.txt
	pyxplot $<

print: $(TXT)
#	$(PIPE) $<

$(TXT) : $(EXE)
	./$< $(DATA) #1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OBJ) $(IMGFL)
