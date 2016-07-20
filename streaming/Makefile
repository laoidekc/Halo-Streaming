CC = mpicc
CFLAGS = -lm

default: halo_streaming
	echo $(CC)

debug: CFLAGS += -g
debug: swan

swan: CC = cc
swan: halo_streaming
	echo $(CC)

archer: CC = cc
archer: halo_streaming
	echo $(CC)

vamp: CC = vtcc
vamp: CFLAGS = -vt:cc mpicc -lm
vamp: halo_streaming
	echo $(CC)

halo_streaming: main.o
	$(CC) $(CFLAGS) -o halo_streaming main.o

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

clean:
	$(RM) halo_streaming *.o *~ streaming_out*.dat streaming_out*.bin exchange_out*.dat exchange_out*.bin halo_streaming.bolt core
