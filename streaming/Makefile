CC = mpicc
CFLAGS = -lm -fastsse

default: halo_streaming

halo_streaming: main.o
	$(CC) $(CFLAGS) -o halo_streaming main.o

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

clean:
	$(RM) halo_streaming *.o *~ out.dat out.bin halo_out.dat halo_out.bin
