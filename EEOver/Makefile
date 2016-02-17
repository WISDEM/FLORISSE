CC      = g++
CFLAGS  = -Wall -g -w 
LDFLAGS = -lgsl -lgslcblas

OBJ = standaloneEEover.o solvers.o zsolve_quartic.o Roots3And4.o

prog: $(OBJ)
	$(CC) $(CFLAGS) -o overlap $(OBJ) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf $(BIN) $(OBJ)
