CC = cc
CFLAGS = -lm

default: extraction_testing

extraction_testing: extraction_testing.o
	$(CC) $(CFLAGS) -o extraction_testing extraction_testing.o

extraction_testing.o: extraction_testing.c
	$(CC) $(CFLAGS) -c extraction_testing.c

clean:
	$(RM) extraction_testing *.o *~ extraction_output*.dat
