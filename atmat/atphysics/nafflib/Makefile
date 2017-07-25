## Makefile for Naff example
#
## Laurent S. Nadolski Synchrotron SOLEIL
#  add with the distribution of M. Gastineau IMCCE

SHELL=/bin/sh

CC = gcc
CCFLAGS=  -O2 

SRC = modnaff.c example.c complexe.c
OBJ  = $(SRC:.c=.o)
INC = modnaff.h complexe.h

SYSLIB = -lm

example : $(OBJ) $(INC)
	@echo "---" $@ "---"
	$(CC)  $(CCFLAGS) $(OBJ) $(SYSLIB) -o  $@

.c.o: $(INC) 
	@echo "---" $@ "---"
	$(CC) $(CCFLAGS) -I. -c $< -o $@

clean: 
	rm -f $(OBJ) $(OBJLIB) example
	
# DO NOT DELETE

