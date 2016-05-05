CC = cc
CFLAGS = -lm

default: injection_testing

injection_testing: injection_testing.o
	$(CC) $(CFLAGS) -o injection_testing injection_testing.o

injection_testing.o: injection_testing.c
	$(CC) $(CFLAGS) -c injection_testing.c

clean:
	$(RM) injection_testing *.o *~ receives_output.dat sends_output.dat
