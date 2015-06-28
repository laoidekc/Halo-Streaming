CC = cc
CFLAGS = -lm

default: halo_streaming
	echo $(CC)

archer: CC = cc
archer: halo_streaming
	echo $(CC)

vamp: CC = vtcc
vamp: CFLAGS = -vt:cc mpicc -lm -fastsse
vamp: halo_streaming
	echo $(CC)

halo_streaming: main.o
	$(CC) $(CFLAGS) -o halo_streaming main.o

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

clean:
	$(RM) halo_streaming *.o *~ out.dat out.bin halo_out.dat halo_out.bin
