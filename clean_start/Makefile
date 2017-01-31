CC = mpicc
CFLAGS = -lm
DEPS = parameters.h
OBJ = main.o parameters.o

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

halo_streaming: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	$(RM) halo_streaming halo_streaming.o* halo_streaming.e* *.o *~ streaming_out*.dat streaming_out*.bin exchange_out*.dat exchange_out*.bin halo_streaming.bolt core
