CC = cc
CFLAGS = -lm

default: basic_streams

basic_streams: basic_streams.o
	$(CC) $(CFLAGS) -o basic_streams basic_streams.o

basic_streams.o: basic_streams.c
	$(CC) $(CFLAGS) -c basic_streams.c

clean:
	$(RM) basic_streams *.o *~ output_basic_streams.dat
