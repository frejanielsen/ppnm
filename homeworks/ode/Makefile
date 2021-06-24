CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat


OUTPUT  =  harm_output.txt SIR_output.txt SIR_output2.txt
EXE   	=  main
OBJ  	=  main.o rung_kutta.o diff_eqs.o
HDRS    =  rung_kutta.h diff_eqs.h
IMGFL   =  harm_plot.png SIR_plot.png SIR_plot2.png
TXT     =  output.txt


.PHONEY: default
default: print $(IMGFL)

print: $(TXT)
	$(PIPE) $<

harm_plot.png: harm_plot.gpi $(HARMOUTPUT)
	pyxplot $<

SIR_plot.png: SIR_plot.gpi $(SIROUTPUT)
	pyxplot $<

SIR_plot2.png: SIR2_plot.gpi $(SIR2OUTPUT)
	pyxplot $<

$(TXT) : $(EXE)
	./$< $(OUTPUT) 1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OUTPUT) $(OBJ) $(IMGFL)
