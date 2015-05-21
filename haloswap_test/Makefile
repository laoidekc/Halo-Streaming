CC = mpicc
CFLAGS = -lm -fastsse

default: halo_swapping

halo_swapping: haloswap.o
	$(CC) $(CFLAGS) -o halo_swapping haloswap.o

haloswap.o: haloswap.c
	$(CC) $(CFLAGS) -c haloswap.c

clean:
	$(RM) halo_swap *.o *~ halo_out.bin halo_out.dat
